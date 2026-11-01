# ECT 共形求解器 2D XZ TE 模支持：项目报告

<!-- SPDX-License-Identifier: BSD-3-Clause -->


| 项 | 内容 |
|---|---|
| 报告日期 | 2026-10-31 |
| 代码基线 | WarpX fork `<repo>`（git 工作区，`master` 尚无提交） |
| 构建目录 | `build_2d`（2D XZ）、`build`（3D） |
| 可执行文件 | `build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB` |
| 相关算例 | `benchmarks/cases/c02_0_pec/` |
| 当前阶段 | **M2a-1b 的 B 侧（2D TE）已完成**：时序缺陷已修复并验证（L2a 逐位 PASS），原 L2 判据前提被证伪并订正为 L2a/L2b（均 PASS）。下一步 L3 + C 侧（E 场对偶共形） |

---

## 0. 摘要

本项目要让 WarpX 的 **ECT（Embedded-boundary Conformal FDTD）** 求解器支持 **2D XZ 的 TE 模**（`Ey` 出平面，`Bx`/`Bz` 入平面）。上游 ECT 在 2D XZ 下只实现 TM 模：`Rhox`/`Rhoz` 被编译期跳过、`face_areas[0,2]` 与 `edge_lengths[1]` 取值为 0、`Bz` 从不更新。

本项目用 **per-unit-y（单位 y 长度）约化** 把 3D 共形面环公式在 2D 下正确退化，使 ECT 框架可直接消费 `Rhox`/`Rhoz`，无需另写 2D 专用离散。代码侧改动分布在 5 个源码文件（几何量约化、`Rho` 放开、`Bz` 放开、稳定化设施、E 更新掩码），均已编译通过；另有两处临时诊断（§4.1 第 13、14 项）。

验证分三层。**L1 几何判据已通过**：网格对齐 PEC 平面下 `Sx={0,dz}`、`Sz={0,dx}`、`ly={0,1}` 精确成立。原 L2 判据（对齐 PEC 上 ECT 必须逐位等于 Yee）**已订正**（见 §14）：

- **L2a（规则域逐位退化）**：EB 面移出域外、无任何 cut cell 时，ECT 与 Yee **201 步逐位一致**（`max|d| = 0`，含 PEC 域壁反射段）→ **PASS**。
- **L2b（对齐壁面位置）**：对齐 PEC 上两者不再逐位相等，差异被完全归因于"贴面 cut cell"的掩膜约定；以反射峰时刻定量：ECT 21.700 ns vs 解析 21.682 ns（误差 0.018 ns），Yee 19.577 ns（等效壁面内移一层 dz）→ **PASS**。

导致原 L2 失败的初始化时序缺陷已修复：`EvolveECTRho` 原先在初始 E 场写入 `Efield_fp` 之前调用，使首个半步 `B^{1/2}` 缺少 `dt·curl(E⁰)`、整体滞后半个时间步；现改为在 `WarpX::AddExternalFields` 末尾（`ApplyEfieldBoundary` 之后）重算，同时覆盖初始化与演化两条路径（§9.2 方案 A 已落地）。

---

## 1. 项目背景

### 1.1 ECT 的三套更新

| 更新 | 文件 | 依赖几何量 |
|---|---|---|
| E | `FiniteDifferenceSolver/EvolveE.cpp` | `face_areas`、`edge_lengths` |
| B | `FiniteDifferenceSolver/EvolveB.cpp` → `EvolveBCartesianECT` | `face_areas`、`area_mod`、`ECTRhofield`、`Venl`、`flag_info_cell` |
| ρ（辅助环量场，**非**电荷密度） | `FiniteDifferenceSolver/EvolveECTRho.cpp` → `EvolveRhoCartesianECT` | `face_areas`、`edge_lengths` |

`ECTRho` 是"E 场沿面周边的积分"，供 B 更新做对偶环量：

```
Bx -= dt·Rhox,  Rhox = ( Ey·ly|_j − Ey·ly|_{j+1} + Ez·lz|_k − Ez·lz|_{k+1} ) / Sx
Bz -= dt·Rhoz,  Rhoz = ( Ex·lx|_j − Ex·lx|_{j+1} + Ey·ly|_{i+1} − Ey·ly|_i ) / Sz
```

`S*` 为面开放面积，`l*` 为边开放长度，均在 EB 初始化时一次性算好。

### 1.2 上游的 2D 限制

1. `EvolveRhoCartesianECT` 中 `Rhox`/`Rhoz` 被 `#ifndef WARPX_DIM_XZ` 排除，XZ 只算 `Rhoy`；
2. XZ 下 `face_areas[0]`、`face_areas[2]` 为 **0**（AMReX 2D 的 `area_frac` 语义是"边长分数"，不是面积）；
3. XZ 下 `edge_lengths[1]` 为 **0**。

三者叠加 ⇒ `Rhox`/`Rhoz` 在 XZ 下是 `0/0`，TE 模无法成立，上游 2D 仅支持 TM。

### 1.3 目标决策

| 决策项 | 结论 |
|---|---|
| 目标维度 | 2D 与 3D 都要生产级 TE |
| 实现路线 | "B（几何约化）+ C（E 侧对偶共形 Ampere）"合并为同一套**对偶几何**设计 |
| 代码边界 | 允许改 fork，不必等上游 |
| 小胞元稳定化 | 扩展/借用策略（`Venl` + `area_mod`） |
| RZ | 不做 |
| AMR | 必须覆盖 |
| 一致性 | 2D TE 须与 3D 逐点一致 |

---

## 2. 实现目标与验收判据

**总目标**：在不改变 3D 数值行为的前提下，使同一套共形面环离散在 2D XZ 下正确退化，让 ECT 在 2D 同时支持 TM 与 TE，并在 cut cell 下保持亚网格精度。

| 层级 | 判据 | 内容 | 状态 |
|---|---|---|---|
| **L1 几何** | 网格对齐时几何严格退化 | `Sx={0,dz}`、`Sz={0,dx}`、`Sy={0,dx·dz}`；`l*={0,1}`；无 face 扩展 | **PASS** |
| **L2a 规则域逐位退化** | 无 cut cell 时场值与 Yee 逐位一致 | EB 面移出域外；`fe`/`fp1`/`fp2` 201 步 | **PASS**（`max\|d\| = 0`） |
| **L2b 对齐壁面位置** | 网格对齐 PEC 上等效壁面落在真实位置 | 反射峰时刻 vs 解析；与 Yee 对比 | **PASS**（+0.018 ns vs Yee 的 −2.1 ns） |
| **L3 斜入射** | 斜面 PEC 反射：振幅 / 角度 / 偏振 | 反射束上 6 个探针（s = 3..8 m）+ 解析到时 + 斜率测角 | **45°/60° × TE/TM 四组 PASS**（平均 `\|r∓1\| = 0.026–0.051`，交叉偏振 0，斜率误差 0.1–2.4%）；**30° 掠入射 FAIL**（斜率 10.5%、TM 横向性 0.45）；见 `CASES_2026-10-31.md` §7quater.2 |
| **L4 圆柱散射** | PEC 圆柱 R=1 的边界/阴影/方位分布 | 8 方位探针 + 柱心 | **PASS**：柱内 0，镜面/阴影 TE 5.35、TM 1.42，前后不对称 0.81/0.30 |

L1 必要不充分；**L2a 才是等价性证明**（规则域上共形格式必须逐位退回 Yee），L2b 用来定量刻画"贴面 cut cell"下两种格式的必然差异（§14.3）。L3 须待反射束探针布置修正后建立。

可复现入口：`benchmarks/cases/c02_0_pec/run_l2.sh`（依次跑 4 个算例并调用 `verify_l2a.py` / `verify_l2b.py`）。

**非目标**：RZ 的 TE；改动 3D 数值路径（3D 必须 bitwise 不变）；本阶段不做性能优化。

---

## 3. 技术方案：per-unit-y 约化

### 3.1 核心洞察

2D XZ 可视为 3D 在 y 方向平移不变的退化，故 2D 几何量应解释为**单位 y 长度**的量：3D 的"面"→ 开放**线段长度**；3D 中沿 y 的"棱"→ **网格节点**，其长度即节点开放度（0/1）。

按此解释，3D 分子中 `Ez` 项因 y-不变性在两个 y 边界严格相消，剩余 `Ey` 项正好给出对 z（或 x）的差商。

### 3.2 映射表

| 几何量 | 上游 2D | per-unit-y 正确值 | 含义 | 网格位置 |
|---|---|---|---|---|
| `face_areas[0]`（x 法向） | 0 | `getAreaFrac()[0] × dz` | x 法向线段开放长度 | 恰为 `Bx(i, j+½)` |
| `face_areas[2]`（z 法向） | 0 | `getAreaFrac()[1] × dx` | z 法向线段开放长度 | 恰为 `Bz(i+½, j)` |
| `edge_lengths[1]`（y 向） | 0 | 节点开放度（流体 1 / 覆盖 0） | 3D y 棱的 2D 退化 | `Ey(i, j)` |

### 3.3 退化验证（解析）

`Rhox = (Ey·ly|_j − Ey·ly|_{j+1} + Ez·lz|_k − Ez·lz|_{k+1}) / Sx`：

- `Ez` 项在 y-不变性下**严格相消**；
- 剩余 `(Ey·ly|_j − Ey·ly|_{j+1})/Sx`，规则区 `ly=1`、`Sx=dz` ⇒ `−(Ey_{j+1}−Ey_j)/dz`；
- 与 Yee 的 `Bx += dt·UpwardDz(Ey)` 一致（符号被 `B -= dt·Rho` 吸收）。

`Rhoz` 同理给 `−(Ey_{i+1}−Ey_i)/dx`。**这保证无 cut cell 时 ECT 在代数层面完全等价于 Yee**，是 L1→L2 的理论依据。

### 3.4 不可省略的配套修改

即使几何量补齐，`EvolveBCartesianECT` 的循环上限是 `idim < AMREX_SPACEDIM`，2D 下只覆盖 `idim = 0,1`，即 **`Bz`（idim=2）永不更新**；须扩为 `idim < 3`，否则 `Rhoz` 算了也没人消费。

---

## 4. 代码改动清单（以源码 diff 为准）

> **修订（2026-10-31 12:45）**：本节按 `diff -u ../WarpX_changed/<file> <file>` 重写。
> 旧版把上游既有的 `MarkUpdateECellsECT`（声明 + 调用）误列为"新增"，并遗漏了
> `WarpXFaceExtensions.cpp` 的稳定化改动；两处均已订正。
>
> 基线：`~/WarpX_changed`（clone `a389dc54a`）；对比：`<repo>`。
> 改动涉及 **5 个源码文件**；行号为 **fork 当前工作树**（`Source/` 下）。

### 4.1 实际改动

| # | 文件 | 函数 / 行号 | 改动（diff 事实） | 状态 |
|---|---|---|---|---|
| 1 | `EmbeddedBoundaryInit.cpp` | `MarkUpdateECellsECT`（:233；hunk 250–280） | Ey 的 XZ 掩膜改由 `ly`（节点开放度）判定，替换上游的 `lx`/`lz` 四边判据 | 完成 |
| 2 | `EmbeddedBoundaryInit.cpp` | `MarkUpdateBCellsECT`（:299；hunk 320–363） | XZ 下 Bx/Bz 的掩膜改用 `face_areas[0]/[2]`（Bx/Bz 所在对偶面），不再用 `lz`/`lx` | 完成 |
| 3 | `EmbeddedBoundaryInit.cpp` | `MarkExtensionCells`（:370；hunk 391–471） | 去掉 XZ 对 `idim=0/2` 的跳过；补 XZ 的 `S_stab`（dim=0 用 `ly(j), ly(j+1), lz`；dim=2 对称）；**并对 `idim=0/2` 直接做 BCK**：`S < S_stab` 时把 `S` 抬到 `S_stab` 并标 `bck_stabilized`，**不进入扩展/借用** | 完成 |
| 4 | `EmbeddedBoundaryInit.cpp` | `ComputeEdgeLengths`（:498；hunk 508–590） | XZ 跳过 `idim=1` 的 `setVal(0)`；在 `lx`/`lz` 之后填 `edge_lengths[1]` = 节点开放度（4 条边任一 covered → 0） | 完成 |
| 5 | `EmbeddedBoundaryInit.cpp` | `ComputeFaceAreas`（:592；hunk 600–659） | XZ 下用 `getAreaFrac()` 填 `face_areas[0]/[2]`（regular→1、covered→0、cut→fraction）；regular/covered box 单独处理（MultiCutFab 无数据） | 完成 |
| 6 | `EmbeddedBoundaryInit.cpp` | `ScaleAreas`（:714；hunk 722–732） | 去掉 XZ 对 `idim=0/2` 的跳过；`full_area = cell_size[(idim+1)%3]*cell_size[(idim+2)%3]`，2D 下 `cell_size[1]=1` → `dz`/`dx` | 完成 |
| 7 | `WarpXFaceExtensions.cpp` | `ComputeSStab`（:84；hunk 96–110、119–133） | 新增 XZ 的 `dim=0/2` 分支：dim=0 用 `0.5·max{ly(j)·dz, ly(j+1)·dz, lz·dy}`，dim=2 对称 | 完成 |
| 8 | `WarpXFaceExtensions.cpp` | `CountExtFaces`（:32；hunk 40–70） | XZ 循环 `idim=1..` → `0..2`；`ReduceIntSum(&sums(0), 3)` | 完成 |
| 9 | `WarpXFaceExtensions.cpp` | `shrink_borrowing`（:240；hunk 258） | `idim < AMREX_SPACEDIM` → `idim < 3` | 完成 |
| 10 | `WarpXFaceInfoBox.H` | :56 | `int vecs_size;` → `int vecs_size = 0;` | 完成 |
| 11 | `EvolveB.cpp` | `EvolveBCartesianECT`（:272） | 分量循环 `idim < AMREX_SPACEDIM` → `idim < 3`（Bz 参与更新） | 完成 |
| 12 | `EvolveECTRho.cpp` | `EvolveRhoCartesianECT`（:113–152） | XZ 下 `Rhox`/`Rhoz` 不再被 `#ifndef` 排除；按 per-unit-y 差分计算 | 完成 |
| 13 | `EvolveB.cpp` | 原 :123–131 | **临时诊断** `DBGB`（每半步全场 `max` 归约） | **已删除** |
| 14 | `EvolveECTRho.cpp` | 原 :159–172 | **临时诊断** `DBGRHO`（每步采样打印） | **已删除** |
| 15 | `WarpXInitData.cpp` | `InitLevelData`（原 :1320–1332） | **删除**过早执行的 `EvolveECTRho`（在初始 E 写入 `Efield_fp` 之前，得到 `ρ ≡ 0`） | 完成 |
| 16 | `WarpXInitData.cpp` | `AddExternalFields`（末尾） | **新增** `EvolveECTRho` 重算（`ApplyEfieldBoundary` 之后），同时覆盖初始化与演化路径、覆盖所有层 | 完成 |

### 4.2 上游既有、**不是** fork 的改动（旧版 §4 误列为"新增"）

| 项 | 位置 | diff 事实 |
|---|---|---|
| `MarkUpdateECellsECT` 声明 | `EmbeddedBoundaryInit.H:75` | **SAME**（与上游逐字节一致） |
| `MarkUpdateECellsECT` 调用 | `WarpXInitData.cpp:1630` | **SAME**（与上游逐字节一致） |
| `EvolveBCartesianECT` 调度 | `EvolveB.cpp:116-118` | 上游既有（此 hunk 无差异） |
| `init_borrowing` 的 0..2 循环 | `WarpXFaceExtensions.cpp:234` | 上游既有 |

> **结论**：fork 只改了 `MarkUpdateECellsECT` 的**函数体**（XZ 分支，见 4.1 第 1 项）；
> 函数声明与调用点均未改动。旧版把它写成"新增"，**多算了 3 项**（旧版第 4、5、6 项），
> 同时**漏掉 6 项**（§4.1 第 2、3、7、8、9、10 项），会误判改动面。

### 4.3 尚未改动 / 设计如此（对应 `M2A_SOURCE_REVIEW_FINDINGS.md` §4）

| 项 | 位置 | 事实 | 性质 |
|---|---|---|---|
| one-way / eight-ways 扩展循环仍跳过 `idim=0/2` | `WarpXFaceExtensions.cpp:671`、`:831` | TE 分量（Bx/Bz）**按设计不做扩展/借用**，小面直接由 §4.1 第 3 项的 BCK 稳定 | **设计如此**，不是遗漏 |
| `ComputeFaceExtensions` 的 BCK 对 0/2 被 `#if !defined(WARPX_DIM_XZ)` 禁用 | `:612`、`:628` | 同上：0/2 的 BCK 已在 `MarkExtensionCells` 内完成，此处不需要 | **设计如此** |
| `EvolveB` 借用索引 `k + vec(1)` 未加 XZ 分支 | `EvolveB.cpp:313/322/345/354` | 当前 `flag_ext_face=0`（0/2 不进入借用），该分支不会执行；**仅当以后为 0/2 启用借用时才会触发** | **潜在**，当前不可达 |
| E 侧共形 Ampere | 无 `EvolveECartesianECT`（`EvolveE.cpp` 仍为掩膜 Yee） | E 侧尚未开始 | **待做**（M2a-1b 之后） |
| `ScaleEdges` 的 XZ `idim=1` 跳过 | `EmbeddedBoundaryInit.cpp:699` | 未改；`cell_size[1]=1` 且 `ly∈{0,1}`，等价于不缩放 | 无影响 |

> 与 `M2A_SOURCE_REVIEW_FINDINGS.md` §4 的差异：该清单把 one-way/eight-ways/BCK 的 XZ 缺口列为"未改"，
> 那是针对"E/B 两侧都做 ECT 式借用"的完整方案；**本 fork 的第一版对 TE 面（Bx/Bz）选择了 BCK-only**
> （`MarkExtensionCells` 内直接抬面积），因此这几项是设计选择而非遗漏。若后续要改为 ECT 借用，
> 需同时改这三处 + `EvolveB` 的借用索引。

### 4.4 说明

- 第 1–6 项集中在 `EmbeddedBoundaryInit.cpp`（几何量 + 掩膜 + 扩展分类），是 "per-unit-y 约化" 的主体；
- 第 7–9 项是"让小胞元稳定化设施在 2D 可用"的配套（旧版 §4 缺失；对应 M2A findings §4 第 6、7 项）：
  `ComputeSStab<0>/<2>` 供第 3 项的 TE 面 BCK 使用；`CountExtFaces` 的 0..2 统计用于警告输出；
  `shrink_borrowing` 扩到 3 是为 `Bz` 的 borrowing 预留（BCK-only 下暂未用到，但保持数组一致）；
- 第 10–12 项是消费端（`vecs_size` 初始化、Bz 循环、Rho）；
- 第 13–14 项是调试期一次性探针，收尾必须移除（§10）。

---

## 5. 工程设施

### 5.1 `benchmarks/cases/c02_0_pec/`

| 文件 | 用途 |
|---|---|
| `inputs_ect_2d_degen` | **L1/L2b**：ECT + 网格对齐 PEC 平面，64×64，CFL 0.4，200 步 |
| `inputs_yee_2d_degen` | 上述的 Yee 参考（仅求解器不同） |
| `compare_degen.py` | 两算例逐分量对比（历史；两者不再逐位相等，见 §14.3） |
| `inputs_ect_2d_noeb` / `inputs_yee_2d_noeb` | **L2a**：EB 面移出域外的规则域对照（零 cut cell） |
| `verify_l2a.py` | L2a 判定：`fe` + `fp1`/`fp2` 逐位一致（默认容差 1e-14） |
| `verify_l2b.py` | L2b 判定：入射/反射峰时刻 vs 解析，ECT 须优于 Yee |
| `run_l2.sh` | L2a + L2b 一键复现（4 个算例 + 两个判定脚本） |
| `inputs_2d` | **L3**：ECT + 30° 斜面，384×384，PML，800 步，入射/反射双探针 |
| `analyze.py` | 斜入射反射系数 |
| `inputs_2d_tm` / `inputs_tm*` / `analyze_tm.py` / `run_tm_sweep.sh` | TM 基线与网格收敛扫描 |
| `A_TM_BASELINE.md` | TM 基线（PEC 反射 `r = −1`） |
| `B_2D_TE_RESEARCH.md` | TE 调研与 per-unit-y 决策依据 |
| `ROOTCAUSE_ECT_2D.md` | 早期"XZ 只算 Rhoy"的根因调查 |
| `dump.py` / `dump_tm.py` | 场数据提取 |

### 5.2 退化算例关键参数

```
max_step = 200 ; amr.n_cell = 64 64 ; amr.max_grid_size = 64  → 单 grid（排除 ghost 交换）
geometry.prob_lo/hi = -8 -8 / 8 8 ; boundary = pec pec（两方向）
warpx.cfl = 0.4  →  dt = 2.358654337e-10 s（含 2D 的 1/√2 因子）
warpx.eb_implicit_function = "z"   →  PEC 半空间，平面与网格严格对齐
Ey = E0·exp(-((z-zc)/w)²)·cos(kz·z),  E0=1e6, zc=-5, w=0.6, kz=2.0944 rad/m,  Bx=-Ey/c
```

**设计要点**：`zc=-5`、`w=0.6` ⇒ `z=0` 处包络 `e^{-69} ≈ 0`，初始条件与导体边界相容；`Ey` 严格与 x 无关，精确解中 `Bz ≡ 0`（PEC 壁平行于 z，不破坏该对称性），因此**任何 `Bz ≠ 0` 都是极敏感的诊断信号**。

---

## 6. 验证进展与调试图时间线

| 阶段 | 实验 | 结果 | 结论 |
|---|---|---|---|
| 1 | L1 几何判据（`S`、`l` 的 min/max/sum） | 仅取 `{0, dz}`、`{0, dx}`、`{0,1}`，无 face 扩展 | **PASS** |
| 2 | 上游回归 `rotated_cube` | 与既有 checksum 一致 | 未破坏上游行为 |
| 3 | 探针逐点对比 | 第 1 步差异已达 1e-9 量级，随后单调放大 | 系统性差异，非噪声 |
| 4 | 移除 EB | ECT 直接 abort（ECT 依赖 EB） | 实验不可行 |
| 5 | **EB 移出域外**（`z - 100`）构造规则域 | 差异依然存在 | **排除 EB/cut cell 因素** |
| 6 | 全场 `FieldEnergy`（步长 1） | 第 0 步全同，**第 1 步即分叉** | 一步内发生，非累积误差 |
| 7 | CFL 扫描（0.4 → 0.1） | 探针相对差 5.8% → 1.4%（比 ≈ 4.1） | 误差**一阶于 dt** ⇒ 半步错位类 |
| 8 | **内核替换**：ECT 临时走 Yee 的 B 内核 | 与 Yee **逐位完全一致**（每步 `rel = 0.000e+00`） | **差异 100% 来自 `EvolveBCartesianECT` 或其输入 `Rho`**；E 更新、掩码、几何均无误 |
| 9 | `DBGRHO`/`DBGB` 幅值打印 | 首个半步前 `Rho ≡ 0`；初始化时 `Ey = 0` | **根因：初始化时序** |

阶段 8 是决定性的：把搜索空间从"整个求解器"收缩到"一个函数 + 一个输入数组"，阶段 9 随即命中。

---

## 7. 根因分析

### 7.1 现象与指纹

规则域、单 grid、无 cut cell，`FieldEnergy` 逐位比较：

| step | 量 | ECT | Yee | 相对差 |
|---|---|---|---|---|
| 0 | total (J) | 4.05301798772508e+01 | 4.05301798772508e+01 | 0 |
| 1 | total | 4.07538790237160e+01 | 4.05290908193397e+01 | **+0.55%** |
| 1 | E | 2.10499360613235e+01 | 2.02297715352119e+01 | **+4.05%** |
| 1 | B | 1.97039429623925e+01 | 2.02993192841279e+01 | **−2.93%** |

指纹：第 0 步全同、第 1 步即分叉；**E 偏大、B 偏小且符号相反** —— 相位错位的典型特征（等相位面上 E/B 相差 π/2），而非幅度或符号错误。

B 场半步步末最大值进一步印证：

| 时刻 | 求解器 | `\|Bx\|max` | `\|Bz\|max` |
|---|---|---|---|
| 第 1 半步末 | ECT | 2.258472009e-3 | **0（严格）** |
| 第 1 半步末 | Yee | 2.287963033e-3 | 3.434207869e-4 |
| 第 2 半步末 | ECT | 2.247439273e-3 | 3.465555913e-4 |
| 第 2 半步末 | Yee | 2.348366e-3 | 6.676228873e-4 |

两条关键证据：

1. **ECT 第 2 半步末的 `Bz`(3.4656e-4) ≈ Yee 第 1 半步末的 `Bz`(3.4342e-4)** ⇒ ECT 比 Yee 慢整整一个半步；
2. ECT 第 1 半步末的 `Bx` 恰为初始条件幅值、`Bz` 严格为 0 ⇒ 该半步**完全没有执行 `curl E` 贡献**。

### 7.2 根因

`EvolveECTRho` 在初始化阶段被调用时 **E 场仍全零**，故 `ECTRhofield ≡ 0`；第一个半步的 B 更新执行的其实是 `B -= dt·0`，B 原地不动。

初始化时序（`Source/Initialization/WarpXInitData.cpp`）：

```
WarpX::InitData()                                   // :797
 ├─ InitFromScratch → MakeNewLevelFromScratch
 │   └─ WarpX::InitLevelData(lev)                   // :1232
 │       ├─ 常量外场 setVal
 │       ├─ InitializeEBGridData(lev)               // :1286  几何量就绪
 │       ├─ [parse 分支] EvolveECTRho(Efield_fp)    // :1324  ← E 仍为 0 ⇒ Rho = 0（缺陷点）
 │       ├─ ComputeExternalFieldOnGridUsingParser(Efield_aux / Efield_cp)   // 仅 lev>0
 │       ├─ [lev>0] EvolveECTRho(Efield_cp)         // :1351  顺序正确
 │       └─ LoadExternalFields(lev)                 // 仅写入 *_external
 └─ for lev: AddExternalFields(lev)                 // :922 → :963
     ├─ Efield_fp += Efield_fp_external             // ← 初始 E 此刻才进入 fp
     └─ ApplyEfieldBoundary(...)
```

即：**`:1324` 的 ρ 初始化发生在初始 E 注入 `Efield_fp` 之前**。`lev > 0` 的 `Efield_cp` 分支（`:1352`）因为其填场在同函数内先执行，反而是正确的 —— 这也解释了为何该缺陷只在 `lev = 0` 上暴露。

### 7.3 数值后果推导

设初值 `(E⁰, B⁰)`，半步 B 更新应为

```
B^{1/2} = B⁰ − (dt/2)·curl(E⁰)
```

ECT 实际得到 `B^{1/2} = B⁰`（因 `Rho = 0`），缺失项量级 `(ω·dt/2)·|B|`：

- `ω = c·kz = 6.2832e8 rad/s`，`dt = 2.3587e-10 s` ⇒ `ω·dt/2 = 0.074`，即**约 7.4% 的 B 幅值**；
- 随后 `E¹ = E⁰ − dt·c²curl(B^{1/2})`，缺失项使 E 偏大、B 偏小，与实测 +4.05% / −2.93% 同量级同符号；
- 半个半步的相位滞后在固定探针上表现为 `τ·dE/dt` 型误差，**正比于 dt**，与 CFL 扫描的 5.8% → 1.4% 完全吻合。

### 7.4 附带发现

`WarpXEvolve.cpp:330` 在"外场随时间变化"时也会在演化中调用 `AddExternalFields`。若按同样逻辑修改 E，则 `ECTRhofield` 同样会变陈旧 —— **修复应放在 `AddExternalFields` 末尾，可一并覆盖初始化与演化两条路径**。

---

## 8. 当前阶段状态矩阵

| 工作项 | 目标 | 状态 | 证据 / 说明 |
|---|---|---|---|
| 几何 per-unit-y 约化 | `S`、`l` 在 2D 有正确取值 | **完成** | `Sx={0,dz}`、`Sz={0,dx}`、`ly={0,1}` 实测 |
| `Rhox`/`Rhoz` 在 XZ 可算 | 去掉编译期屏蔽 | **完成** | 探针显示 `Rhox=1.42e6` 等非零值 |
| `Bz` 参与 ECT 更新 | `idim < 3` | **完成** | `DBGB`（已删除）曾显示 `Bz` 非零 |
| E 更新掩码 `MarkUpdateECellsECT` | 上游既有；fork 仅改 XZ 判据（`ly` 节点开放度） | **完成（XZ 语义）** | 函数/调用为上游既有（§4.2）；内核替换实验证明 E 侧无差异 |
| 小胞元稳定化设施（`ComputeSStab`/`CountExtFaces`/`shrink_borrowing`） | 2D 可用 | **完成** | §4.1 第 7–9 项；TE 面按设计只做 BCK（§4.3），不做借用 |
| L1 几何判据 | 严格退化 | **PASS** | §6 阶段 1 |
| 根因定位 | 找出唯一差异源 | **完成** | §7.2，内核替换实验逐位一致 |
| 修复实施（方案 A） | 时序修正 | **完成** | `WarpXInitData.cpp`：删 `:1324` 旧调用、`AddExternalFields` 末尾重算（§4.1 第 15–16 项） |
| L2a 规则域逐位退化 | `max\|d\| = 0` | **PASS** | 201 步 `fe`/`fp1`/`fp2` 全为 0（§14.2） |
| L2b 对齐壁面位置 | 等效壁面在真实位置 | **PASS** | ECT 反射峰 21.700 ns vs 解析 21.682 ns；Yee 19.577 ns（§14.3） |
| L3 斜入射判据 | `r → −1` | **未建立** | 见 §12.1 |
| 临时诊断清理 | 移除 `DBGB`/`DBGRHO` | **完成** | `grep DBGB\|DBGRHO Source/` 无命中 |
| 上游回归 | 2D/3D `rotated_cube` 不回归 | **PASS** | 2D 相对差 6.8e-10、3D ≤ 1e-13（§14.4） |

**一句话概括现在阶段**：2D TE 的实现 + 时序缺陷修复均已完成，L2a（逐位退化）与 L2b（壁面位置）双双通过；**B 侧（Faraday/环量）到此结束**，下一步是 L3 判据重建与 C 侧（E 场对偶共形 Ampere）。

---

## 9. 修复方案

### 9.1 方案对比

| 方案 | 做法 | 优点 | 缺点 |
|---|---|---|---|
| **A（推荐）** | 把 `EvolveECTRho` 的初始化调用移到 `WarpX::AddExternalFields` 末尾（`ApplyEfieldBoundary` 之后） | 一次覆盖初始化与演化中改场两条路径；`lev=0` 与 `lev>0` 统一 | 需在 `AddExternalFields` 内取 `edge_lengths`/`face_areas`，多几行样板 |
| B | 保留 `:1324`，在 `InitData` 中 `AddExternalFields` 之后补一次调用 | 改动最小 | 演化路径（`WarpXEvolve.cpp:330`）仍会陈旧 |
| C | 每个 B 半步前重算 ρ | 绝对安全，无需考虑时序 | 每步多一次全网格 kernel，性能损失 |

### 9.2 方案 A 代码骨架

```cpp
// Source/Initialization/WarpXInitData.cpp, WarpX::AddExternalFields 末尾
    if (m_p_ext_field_params->E_ext_grid_type != ExternalFieldType::default_zero) {
        ApplyEfieldBoundary(lev, PatchType::fine, 0.0);

#ifdef AMREX_USE_EB
        // Recompute the ECT circulation field so that it is consistent with the
        // E field that the solver will actually see. The ECTRho initialization
        // inside InitLevelData runs before AddExternalFields, i.e. before the
        // initial E field reaches Efield_fp, and would otherwise stay zero and
        // make the first half-step of the B update a no-op.
        if (EB::enabled() &&
            WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT) {
            m_fdtd_solver_fp[lev]->EvolveECTRho(
                m_fields.get_alldirs(FieldType::Efield_fp, lev),
                m_fields.get_alldirs(FieldType::edge_lengths, lev),
                m_fields.get_mr_levels_alldirs(FieldType::face_areas, max_level)[lev],
                m_fields.get_alldirs(FieldType::ECTRhofield, lev),
                lev);
        }
#endif
    }
```

同时**删除** `InitLevelData` 中 `:1320-1332` 的旧调用（保留 `lev>0` 的 `Efield_cp` 分支，或一并迁移，视 AMR 一致性测试结果而定）。

### 9.3 验收方式（已于本次迭代执行，见 §14）

1. **L2a**：`inputs_ect_2d_noeb` vs `inputs_yee_2d_noeb`（EB 面移出域外），要求 `FieldEnergy` 与 `fp1`/`fp2` 每步 `rel = 0`；
2. **L2b**：`inputs_ect_2d_degen` vs `inputs_yee_2d_degen`（对齐 PEC），要求 ECT 反射峰与解析时刻之差 ≤ 0.5 ns 且优于 Yee；
3. 跑 2D/3D `rotated_cube` 回归，确认 checksum 不变；
4. 再进入 L3 斜入射（探针布置需先修正，§12.1）与 AMR 测试。

> 说明：原判据"对齐 PEC 上 ECT 必须逐位等于 Yee"**不可达**（§14.3），不要再用 `compare_degen.py` 的 bitwise 阈值作为门禁。

---

## 10. 收尾清单

- [x] 删除 `EvolveB.cpp` 的 `DBGB` 打印块
- [x] 删除 `EvolveECTRho.cpp` 的 `DBGRHO` 打印块（含 `i0`/`j0` 采样变量）
- [x] 实施方案 A 并删除 `WarpXInitData.cpp` 的旧调用
- [x] 更新 `inputs_ect_2d_degen` 头部注释（改为解释"贴面 cut cell"与 L2a/L2b 分工）
- [x] 补充 L2 的自动化断言（`verify_l2a.py` / `verify_l2b.py` + `run_l2.sh`）
- [x] 上游回归复算（2D `rotated_cube` 6.8e-10、3D `rotated_cube` ≤1e-13）
- [ ] `pre-commit run -a`（clang-format + Ruff）
- [ ] 清理算例目录下的临时产物：`Backtrace.0.0`、`diags/`、`warpx_used_inputs`

---

## 11. 复现手册

### 11.1 构建

```bash
cd <repo>/build_2d
cmake --build . --target app_2d -j 20
# 3D: cd ../build && cmake --build . --target app_3d -j 20
```

### 11.2 L2a / L2b 判据

```bash
cd <repo>/benchmarks/cases/c02_0_pec
./run_l2.sh              # 4 个算例 + 两个判定脚本，输出到 /tmp/l2
# 或分步：
#   python3 verify_l2a.py /tmp/l2/ect_noeb /tmp/l2/yee_noeb
#   python3 verify_l2b.py /tmp/l2/ect_pec  /tmp/l2/yee_pec
```

`WARPX_BIN` / `RUNS` / `OMP_NUM_THREADS` 均可用环境变量覆盖，默认
`<repo>/build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB`、`/tmp/l2`、`8`。

参考结果（2026-10-31，64×64，CFL 0.4，200 步）：

```
L2a  RESULT: PASS   fe/fp1/fp2 全部 max|d| = 0.000e+00
L2b  ECT: t_ref = 21.700 ns (解析 21.682 ns, err +0.018 ns), r = -0.9188
     YEE: t_ref = 19.577 ns (err -2.105 ns),                    r = -0.8505
     RESULT: PASS
```

### 11.3 全场能量诊断（定位用）

在算例末尾追加：

```
warpx.reduced_diags_names = fe
fe.type = FieldEnergy
fe.intervals = 1
```

### 11.4 决定性定位实验（内核替换）

把 `EvolveB.cpp` 的 ECT 分支临时替换为

```cpp
EvolveBCartesian <CartesianYeeAlgorithm> ( Bfield, Efield, Gfield, lev, dt );
```

若此时 ECT 与 Yee **逐位一致**，即可确认 E 更新与几何量均无误，差异只在 B 内核或其输入 ρ。**该实验完成后必须还原。**

---

## 12. 遗留问题与风险

### 12.1 L3 判据（已重建，2026-10-31）

旧 `inputs_2d` 的"反射探针" `(0,+4)` **在导体内部**，`r = −0.531` 无意义。

已新建 `inputs_2d_l3`：按镜面几何把 6 个探针放在反射束上（`k̂_r = (−sin2α, cos2α)`，
`s = 3..8 m`），并给出解析到达时刻。结果（384²，ECT）：

- 全部探针 `r < 0`（相位反转），导体内部 `max|Ey| = 0`；
- 沿束平均 `r = −0.9876`（`|r+1| = 0.0124`），单点 −0.916..−1.114；
- 384² 与 768² 几乎相同（差值 < 0.0015）⇒ 起伏来自有限斜面/外边界的衍射，非离散误差；
- Yee（PEC 边界，因 PML 与投影清洁器不兼容）起伏形态相同、平均 −1.076，说明该调制非 ECT 特有。

判据与数据见 `L3_REFLECTION_CRITERION.md`。**剩余工作**：加载波后做 2D FFT 谱分解，
直接核对反射角 `2α` 与幅值（当前无载波脉冲无法在 k 空间分离）。

### 12.2 其他风险

| 项 | 说明 |
|---|---|
| AMR | 当前所有验证都在 `max_level = 0`；`face_areas` 的 per-unit-y 约化在 coarse/fine 界面（`area_mod`、`Venl`）尚未验证 |
| 3D 不回归 | 已实测：2D `rotated_cube` 相对差 6.8e-10、3D `rotated_cube` ≤ 1e-13（§14.4）。注意 3D 的 ECT **同样存在**被修复的初始化时序缺陷（初场 E≠0 时才会显现），本次一并修正 |
| 小胞元稳定性 | cut cell 出现极小 `S` 时，`Venl`/`area_mod` 的 BCK-only 策略尚未做定量稳定性扫描（`r = −0.9188` 提示薄包络脉冲的反射幅值仍有 ~8% 偏差，需 L3 分辨率扫描确认是格式耗散还是稳定化损失） |
| 2D/3D 一致性 | 决策要求"2D TE 与 3D 逐点一致"，该测试尚未搭建（需同一物理设置的三维算例） |
| 贴面 cut cell 语义 | ECT（`ly`/`face_areas` 由几何开放度给）与 Yee（`isRegular()`）在"volfrac=1 但被标记为 cut"的单元上判定相反；这是 L2b 的全部来源，**后续任何与 Yee 的对比都要考虑这一层差异**（§14.3） |
| BCK 与 `S_mod` | TE 面按设计只做 BCK：`S < S_stab` 时把 `S_mod` 抬到 `S_stab` 并标 `bck_stabilized`，B 更新走 `B -= dt·Rho`（不带 `S_mod`）。若后续改为借用式扩展，需重审此分支与 `Venl` 记账 |

---

## 13. 附录

### 13.1 术语

| 术语 | 含义 |
|---|---|
| ECT | Embedded-boundary Conformal FDTD，基于 EB 几何的共形有限差分求解器 |
| TE / TM（2D XZ） | TE：`Ey` 出平面；TM：`By` 出平面 |
| `ECTRhofield` | E 场沿面周边的环量，供 B 更新使用，**非**电荷密度 |
| per-unit-y | 把 2D 几何量解释为单位 y 长度的量，使 3D 公式正确退化 |
| `face_areas` | 面开放面积（2D 下为开放线段长度） |
| `edge_lengths` | 边开放长度（2D 下 y 向退化为节点开放度） |
| `flag_info_cell` | 小胞元稳定化所需的 cell 邻域标记 |

### 13.2 关键文件索引

```
Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp      per-unit-y 约化、MarkUpdateECellsECT
Source/EmbeddedBoundary/EmbeddedBoundaryInit.H        MarkUpdateECellsECT 声明（上游既有，未改）
Source/EmbeddedBoundary/WarpXFaceExtensions.cpp        2D S_stab / CountExtFaces / shrink_borrowing
Source/EmbeddedBoundary/WarpXFaceInfoBox.H             vecs_size 初始化
Source/FieldSolver/FiniteDifferenceSolver/EvolveB.cpp          EvolveBCartesianECT、idim<3
Source/FieldSolver/FiniteDifferenceSolver/EvolveECTRho.cpp     EvolveRhoCartesianECT
Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp          E 更新（共用路径）
Source/Initialization/WarpXInitData.cpp                        InitLevelData、AddExternalFields
Source/Evolve/WarpXEvolve.cpp                                  EvolveB 两次调用、AddExternalFields(:330)
benchmarks/cases/c02_0_pec/                                    全部验证算例与历史报告
```

### 13.3 相关历史文档

> **权威副本在 `Docs/research_notes/`**；`benchmarks/cases/c02_0_pec/` 下的同名 md 是算例目录的历史副本，可能滞后。

- `Docs/research_notes/CASES_2026-10-31.md` —— **本轮所有已跑算例的物理情景与结果**（C1–C11，含失败/作废实验）
- `Docs/research_notes/L3_REFLECTION_CRITERION.md` —— L3 判据重建（反射束探针、解析时刻、30°/60°、15° 限制）
- `Docs/research_notes/M2A_CODE_AUDIT.md` —— 第三方代码审计（§0 为审计后的处置状态与新算例）
- `Docs/research_notes/C_SIDE_DESIGN.md` —— C 侧（E 场对偶共形 Ampere）实现清单 + AMR/2D-3D 对照方案
- `Docs/research_notes/ROOTCAUSE_ECT_2D.md` —— 早期"XZ 只算 Rhoy"的根因调查
- `Docs/research_notes/B_2D_TE_RESEARCH.md` —— TE 调研与 per-unit-y 决策（附录为决策后修正）
- `Docs/research_notes/A_TM_BASELINE.md` —— TM 基线（PEC 反射 `r = −1`，含修复后复算）
- `Docs/research_notes/cutcell-research.md` —— M2a 设计规范（§8，含 L2 判据订正注记）
- `Docs/research_notes/M2A_SOURCE_REVIEW_FINDINGS.md` —— 源码级核查清单（D1–D6 待决问题的出处）

---

## 14. 二次迭代记录（2026-10-31 · 修复落地 + L2 判据订正）

> 本节记录 §9 方案 A 落地后发生的全部事实，是 §2 / §8 / §10 更新的依据。

### 14.1 落地的代码改动（4 处）

| # | 文件 | 改动 |
|---|---|---|
| 1 | `EvolveB.cpp` | 删除 `DBGB` 打印块（原来每半步做三次全场 `max()` 归约，严重拖慢且污染日志） |
| 2 | `EvolveECTRho.cpp` | 删除 `DBGRHO` 打印块（含 `i0`/`j0` 采样变量） |
| 3 | `WarpXInitData.cpp` `InitLevelData` | 删除 `:1320–1332` 的 `EvolveECTRho`（早于初始 E 注入 `Efield_fp`，必然得到 `ρ ≡ 0`） |
| 4 | `WarpXInitData.cpp` `AddExternalFields` | 末尾新增 `EvolveECTRho` 重算（`ApplyEfieldBoundary` 之后），条件：`EB::enabled() && solver==ECT && m_fields.has_vector(ECTRhofield, lev)` |

第 4 项的注释说明了必要性：B 更新消费的是 `ECTRhofield` 而不是现算的 `curl E`，一旦环量陈旧，第一个半步就是空操作，整体错位半拍。放在 `AddExternalFields` 里可同时覆盖初始化（`InitData` → `:922`）与演化中改场（`WarpXEvolve.cpp:330`）两条路径、并覆盖所有层；`lev>0` 的 `Efield_cp` 分支（`:1352` 附近）保持不变，因为它的填场就在同函数内先行完成。

### 14.2 L2a 结果：规则域逐位退化 —— PASS

EB 面移出域外（`warpx.eb_implicit_function = "z - 100"`，域为 [−8,8]²）⇒ 全域无 cut cell。ECT vs Yee，201 步：

```
--- fp1 (201 samples) ---            --- FieldEnergy (201 samples) ---
  OK Ey  peak=6.235032e+05 max|d|=0   OK total peak=4.060209e+01 max|d|=0
  OK Bx  peak=1.944326e-03 max|d|=0   OK E     peak=3.069605e+01 max|d|=0
  (其余分量同)                        OK B     peak=2.761955e+01 max|d|=0
RESULT: PASS
```

要点：

1. **`max|d|` 严格为 0**，不是"舍入量级"——说明在规则域上两者的算术路径也一致，per-unit-y 约化在代数层面完全等价于 Yee；此前"0.55% 偏差"确认全部来自时序缺陷。
2. 该算例中脉冲在 t≈45 ns 反射于 `z=+8` 的 PEC 域壁，**边界路径也被覆盖**（域壁用标准 `ApplyEfieldBoundary`，两种求解器共用）。
3. 由此可把 bulk 与 EB 界面解耦：剩余的 ECT/Yee 差异**只可能**来自 EB 掩膜/贴面几何。

### 14.3 原 L2 判据不可达：贴面 cut cell 的掩膜约定不同

**现象**：对齐 PEC 上 ECT 与 Yee 从第 ~28 步开始分叉（`FieldEnergy` 相对差 2e-7 → 4e-3，场差最大 0.35 @ z≈−0.4），且总能量始终几乎相同、E/B 配比差异很大 —— 是**相位/共振频率**差异，不是能量不守恒。

**证据链**：

| # | 实验 | 结果 |
|---|---|---|
| 1 | 打印两种求解器的 `m_eb_update_*` 掩膜 | ECT：`Ey` 节点 `j=0..31` 更新、`j≥32` 冻结；Yee：只到 `j=0..30` ⇒ **Yee 多冻结一层** |
| 2 | 打印 EB 单元标记 | `j=31`（z∈[−0.25,0]，贴面单元）`volfrac=1`、`areafrac=1`，但 `isRegular()=false` |
| 3 | `MarkUpdateCellsStairCase` 源码 | 判据是 `!flag.isRegular()` ⇒ 把贴面单元当封闭 ⇒ 冻结 `j=31` 及其 Bx 行 |
| 4 | ECT 几何 | `ly`/`face_areas` 由 `edge_cent`/`areaFrac` 给 ⇒ 贴面单元按**开放**处理（`Sx=dz`、`ly=1`） |
| 5 | 反射峰时刻（探针 z=−1.5，解析：壁在 z=0 时 21.682 ns，壁内移一层时 20.014 ns） | **ECT 21.700 ns（+0.018 ns）**，Yee 19.577 ns（−2.105 ns） |
| 6 | 反射幅值比 | ECT `r = −0.9188`，Yee `r = −0.8505` |

**结论**：

- 对齐面**不产生"无 cut cell"的规则几何**：AMReX 会把贴面单元标成 cut（volfrac 仍为 1）。原判据"对齐 PEC 上必然逐位等于 Yee"的前提不成立。
- 两者等效壁面相差一层 `dz`：ECT 落在真实位置 `z=0`，Yee 落在 `z=−dz`。差异在物理上是 **O(dx) 的阶梯近似误差**，不是 fork 的缺陷。
- 因此 L2 被拆成 L2a（规则域逐位，证明格式等价）+ L2b（壁面位置，证明 ECT 更贴真几何）。原 `compare_degen.py` 的 bitwise 门禁**作废**（保留作历史工具）。

### 14.4 回归与基线复算

| 算例 | 结果 |
|---|---|
| 2D `embedded_boundary_rotated_cube`（上游 ECT 回归，TM） | 相对差 6.8e-10（Ex/Ez/By 全部吻合，Bx/Bz/Ey 仍严格 0）→ 2D TM 未受影响 |
| 3D `embedded_boundary_rotated_cube` | 相对差 ≤ 1e-13 → 3D 路径未受影响（该算例初场 E≡0，修复对它是恒等操作） |
| 2D TM 基线 `inputs_2d_tm_flat`（n=384） | `r = −0.999361`（修复前 −0.999295），入射 10.024 ns / 反射 43.488 ns 不变 → 修复对 TM 只有二阶影响，方向为更好 |

### 14.5 新增的验证资产

| 文件 | 作用 |
|---|---|
| `inputs_ect_2d_noeb` / `inputs_yee_2d_noeb` | L2a 算例对（EB 移出域外 + `fp1`/`fp2`/`fe`） |
| `verify_l2a.py` | L2a 自动判定（默认容差 1e-14，实测 0） |
| `verify_l2b.py` | L2b 自动判定（入射/反射峰、三条断言） |
| `run_l2.sh` | 一键复现 4 个算例 + 两个判定 |

### 14.6 下一步

1. ~~**L3 判据重建**~~ —— **已完成**：`inputs_2d_l3` 沿反射束 6 探针，ECT 沿束平均 `r = −0.9876`（`|r+1| = 0.0124`），384²/768² 一致；见 `L3_REFLECTION_CRITERION.md`。
2. **C 侧**（`EvolveECartesianECT` + `dual_face_area`/`dual_edge_length` + 对偶 `S_stab`）：实现清单、2D 逐分量公式、退化判据已写入 `C_SIDE_DESIGN.md` §4/§9，**代码尚未开工**。
3. **2D/3D 一致性**：对照设置（1 层 y + 周期 + `dy=1` + 同 dx/dz/dt + 同探针）与容差讨论见 `C_SIDE_DESIGN.md` §7。
4. **AMR 多层几何**：改动点（per-level 初始化、扩展链、粗-细一致化）见 `C_SIDE_DESIGN.md` §6。
5. **审计反馈处置**（`M2A_CODE_AUDIT.md` §6 → §0）：
   - **P1 `verify_l2b.py` 判据** 已修：`dz`/`zc` 自动解析、改用往返时间、`r` 降级为信息项；
     64²/128²/256² 全部 PASS（新增 `run_l2b_scan.sh`）。
   - **P1 C 侧**：按用户指示暂缓。
   - 新增算例：45° 稳定性（`inputs_ect_2d_45deg`）、3D TE（`inputs_3d_te_pulse`）、
     60° 斜入射（`inputs_2d_l3_60deg`）与角度生成器 `make_l3_angle.py`。
   - **未改**：TE 面仍 BCK-only、`S` 就地覆盖、`EvolveB` 借用索引（不可达）、初始化修复覆盖面
     （仅 external-E 路径）——均作为已知限制记录，待 C 侧或后续里程碑处理。
6. **多偏振 × 多角度 + 圆柱 + 长跑**（`CASES_2026-10-31.md` §7quater）：
   - 45°/60° × TE/TM 四组：振幅（2.6–5.1%）、偏振（交叉分量 0）、角度（斜率 0.1–2.4%）全部通过；
   - 30° 掠入射不达标（斜率 10.5%，TM 横向性 0.446）⇒ 当前精度不足以支撑掠入射定量判据；
   - PEC 圆柱：柱内 0、阴影/镜面比 TE 5.35 / TM 1.42、前后不对称 0.81/0.30，全部通过；
   - 10000 步长跑：TE 能量漂移 3.0e-4，**TM 有 1.9% 能量增长**（未发散，但与 TE 差两个量级）
     —— 指向 `Rhoy` 与 BCK 抬面积的交互，列为 C 侧要解决的首个精度问题。

---

*本报告基于 2026-10-31 的工作区状态编写；§2 / §8 / §10 已按 §14 的实测结果更新。*



