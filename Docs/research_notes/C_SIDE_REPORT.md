# C 侧实施报告：改动位置、效果、残留问题

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 范围：从"开工吧，动手"这条指令之后产生的**全部**改动（提交 `6225846` → `f1c6190`）。
> 分支：`feat/ect-conformal-E`（开发）→ `master`（已合并），tag `c-side-2d-v1`。
> 行号均指当前 `master` 工作树。
> 前置文档：`C_SIDE_DESIGN2.md`（设计）、`C_SIDE_RESEARCH2.md`（开工条件）、`Docs/C_SIDE_REVIEW.md`（审查）、
> `C_SIDE_IMPL_LOG.md`（过程记录，本文是它的汇总版并补齐残留问题清单）。

---

## 1. 起点与终点

| | 开工前（`6225846`） | 现在（`f1c6190`） |
|---|---|---|
| E 更新 | 掩膜 Yee，`EvolveE.cpp:92` 把 `Yee` 与 `ECT` 合并到同一模板 | 新增 `EvolveECartesianECT`，开关打开时走共形内核 |
| 对偶几何 | 不存在 | `A_y`、`L_x`、`L_z` 三套 MultiFab，初始化时由 `volfrac` 组装 |
| 开关 | 无 | `algo.ect_conformal_E = 0/1/2` |
| 米氏爆炸 | 归因于"E 侧阶梯"，阻塞 | **归因于多 grid**，单 grid 可用；收敛阶基线已测出 |
| 版本控制 | 只有文档提交 | 代码 + 算例 + 文档全部提交，合并到 `master` 并打 tag |

代码规模：**8 个源文件，+256 行**；算例脚本 3 个；文档 5 份。

---

## 2. 改动位置（逐文件）

### 2.1 数据结构与开关

| 文件 | 位置 | 内容 |
|---|---|---|
| `Source/Fields.H` | `:73-75`（枚举，`area_mod` 之后） | `dual_A_y`（节点对偶面积）、`dual_L_x`（沿 x 开放长度，存于 Bx 处）、`dual_L_z`（沿 z，存于 Bz 处） |
| `Source/WarpX.H` | `:202` | `static inline int ect_conformal_E = 0;`（0 关 / 1 常数几何 / 2 真实几何） |

### 2.2 参数解析与内存分配

| 文件 | 位置 | 内容 |
|---|---|---|
| `Source/WarpX.cpp` | `:582` | `pp_algo.query("ect_conformal_E", ect_conformal_E)` |
| 同上 | `:586-596` | 三条断言：必须 `maxwell_solver = ect`；取值只能是 1/2；非 XZ 构建直接 abort |
| 同上 | `:2714-2727`（`AllocLevelMFs` 的 ECT 段） | 按 `Ey_nodal_flag` / `Bx_nodal_flag` / `Bz_nodal_flag` 分配三套，`ngrow = ng_FieldSolver`，初值 1.0（避免未初始化时除零） |

### 2.3 对偶几何计算

| 文件 | 位置 | 内容 |
|---|---|---|
| `Source/EmbeddedBoundary/EmbeddedBoundaryInit.H` | `:152` | `ComputeDualGeometry(dual_A_y, dual_L_x, dual_L_z, eb_fact, cell_size, constant_geometry)` 声明 |
| `Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp` | `:746` | 实现。三分支：`regular`（填 `dx*dz`/`dx`/`dz`）、`covered`（填 0）、`mixed`（读 `volfrac` 组装） |
| 同上 | `:798-810` | 方案 B 公式：`Ay = ¼dxdz·Σv`（四象限）、`Lx = ½dx(v(i-1,j)+v(i,j))`、`Lz = ½dz(v(i,j-1)+v(i,j))` |
| `Source/Initialization/WarpXInitData.cpp` | `:1639-1655` | `InitializeEBGridData` 的 ECT 分支调用；`verbose` 时打印三者的 min/max |

### 2.4 内核

| 文件 | 位置 | 内容 |
|---|---|---|
| `Source/FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H` | `:537` | `EvolveECartesianECT` 声明 |
| `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp` | `:92-96` | 分派：`ECT && ect_conformal_E > 0` → 新内核，放在原有 `Yee || ECT` 分支**之前** |
| 同上 | `:241-327` | 内核实现（2D XZ，非 XZ 走 `:323` 的 abort）。三条更新见 §3.1 |
| 同上 | `:299-318` | covered 保护：`A <= 0` 时把该分量写 0（不是仅跳过），避免脏值残留 |

### 2.5 算例与后处理

| 文件 | 位置 | 内容 |
|---|---|---|
| `benchmarks/cases/c02_0_pec/verify_l2a.py` | `:47` | 加噪声地板：`noise = field_peak * 1e-13`，零分量（TE 的 `Bz`，峰值 1e-19）不再假 FAIL |
| `benchmarks/cases/c02_0_pec/make_mie_cases.py` | `:52` | `amr.max_grid_size = {n}`（原为固定 64），并注明原因 |
| 同上 | `inputs_mie_{ey,by}_{40,80,160,320,640}` | 重新生成（每档 2 行改动） |
| `benchmarks/cases/c02_0_pec/analyze_mie_convergence.py` | `:80-95` | 新增 `pf_time()`，DFT 相位改用 plotfile Header 的物理时间（原用 `j*T/8`） |
| 同上 | `:125-130` | 数值场除以 `my_constants.E0`（原先未归一化，L2 恒为 1e6） |

### 2.6 文档

| 文件 | 内容 |
|---|---|
| `C_SIDE_RESEARCH.md` / `C_SIDE_DESIGN2.md` / `Docs/C_SIDE_REVIEW.md` / `C_SIDE_RESEARCH2.md` | 开工前已存在，本轮只读不改 |
| `C_SIDE_IMPL_LOG.md` | **新增**：过程记录（代码清单、L2a/L3、根因对照实验、收敛阶基线） |
| `HANDOFF.md` | 改写 §3（米氏归因）、索引新增 3 份文档 |

---

## 3. 效果

### 3.1 内核公式（实现即文档）

```
Ex(i,j) += −c²dt [By(i,j) − By(i,j-1)] / Lz(i,j)
Ez(i,j) += +c²dt [By(i,j) − By(i-1,j)] / Lx(i,j)
Ey(i,j) += +c²dt [Bx(i,j)Lx(i,j) − Bx(i,j-1)Lx(i,j-1)
                  −Bz(i,j)Lz(i,j) + Bz(i-1,j)Lz(i-1,j)] / Ay(i,j)
```

规则域（`Lx=dx, Lz=dz, Ay=dx·dz`）逐项退化为 `EvolveE.cpp:173-199` 的 Yee 模板。

### 3.2 L2a（规则域，EB 在域外）

| 配置 | Ey rel | 能量 rel | 结论 |
|---|---|---|---|
| 常数几何（`=1`） | 1.6e-15 | 2.6e-15 | PASS |
| 真实几何（`=2`） | 1.9e-15 | 2.6e-15 | PASS |

几何量自检：`A_y = 0.0625 = dx·dz`、`L_x = L_z = 0.25`，规则域给出**精确常数**。

### 3.3 L3 45° TE（真实几何）

| 判据 | 旧（掩膜 Yee） | 新（共形 E） |
|---|---|---|
| 斜率误差（角度） | −2.39% | **−1.51%** |
| 平均 `\|r+1\|` | 0.0427 | 0.0459 |
| 交叉偏振 | 0 | 0 |

角度判据改善，幅值略差，净效果不显著（一阶几何与阶梯同阶）。

### 3.4 米氏（TE，单 grid，后处理修正后）

| 配置 | 近场 L2 @80² | @160² | 近场阶数 | RCS 阶数 |
|---|---|---|---|---|
| 掩膜 Yee（`=0`） | 0.794 | 0.327 | **1.28** | **1.76** |
| C 侧（`=2`） | 0.932 | 0.440 | 1.06 | 1.77 |

这是项目第一次拿到**可信的**收敛阶基线（此前 L2 恒为 1e6、相位乱跳）。

### 3.5 副作用：两个此前被掩盖的真相

1. **"逐位一致"不可能**：`(Bx−Bx)/dz` 与 `(Bx·dx−Bx·dx)/(dx·dz)` 数学等价但浮点顺序不同，
   实测 rel 1e-15。C-1 的判据应从 `max|d| = 0` 改为 `rel ≤ 1e-14`。
2. **米氏爆炸与 C 侧无关**：`max_grid_size = 64 < n_cell = 80` 导致 4 个 grid，CW 激光源跨 grid 不稳。
   80² 单 grid 下，掩膜 Yee / C 侧 / 带 EB 圆柱全部稳定在 1.2–1.9e6。历史现象（40²/64² 稳定、
   80²/160² 爆炸）正好对应单/多 grid 的分界。

---

## 4. C 侧残留问题

### 4.1 P1（影响结果可信度或可用范围）

| # | 问题 | 现状 | 影响 | 出路 |
|---|---|---|---|---|
| R1 | **掩膜语义未改** | `MarkUpdateECellsECT` 仍用 `ly == 0` 冻结；`ly == 0` 涵盖"三开放一覆盖"的节点，而这些节点的 `A_y > 0` | 部分 cut 节点被冻结，共形几何在这些点上**根本没生效**，可能是 L3/米氏改善不显著的原因之一 | 按 DESIGN2 §2.3.3 把 `Ey` 判据统一为 `A_y <= 0`（S6' 未做） |
| R2 | **无对偶稳定化** | 实测 `A_y` 最小到 0、`L` 最小到 0；cut 点的小值无任何下限保护 | 极小 `A_y` 处 `1/A_y` 增益无界；当前算例未因此爆炸（多 grid 问题掩盖），但曲面加密后风险高 | 实现 `A_stab`（只抬分母），做 `volfrac` 扫描标定 |
| R3 | **几何是一阶** | 方案 B（`volfrac` 双线性），非 level-set 精确积分 | 收敛阶与阶梯同阶（1.06 vs 1.28），C 侧的价值未兑现 | 方案 C：初始化时用 level set 求交积分（`C_SIDE_DESIGN2.md` §1.5），成本只在初始化 |
| R4 | **3D 未实现** | 内核 `#if defined(WARPX_DIM_XZ)`；3D 构建在 `WarpX.cpp:594` 直接 abort | 无法用于 3D（项目最终要求 2D/3D） | 按 DESIGN2 §5 的三分量四象限推广，先 `A` 后 `L` |
| R5 | **AMR 未支持** | `InitializeEBGridData` 只在 `lev == maxLevel()` 初始化；三套数组的分配初值是 `1.0_rt`，低层不会被 `ComputeDualGeometry` 覆盖 | 多层时低层保持 **1.0**，内核 `÷1.0` ⇒ **量纲错、静默错误更新**（不冻结、不 abort、不除零）。**原报告写的"几何为 0 → 除零 → 冻结"是错的**，见 `C_SIDE_REPORT_CHECK.md` §3 | per-level 初始化；在此之前加 `maxLevel()==0` 断言（`C_SIDE_FIX_PLAN.md` F2） |

### 4.2 P2（生产化缺口，当前测试不触发）

| # | 问题 | 现状 | 影响 |
|---|---|---|---|
| R6 | **PML 区未排除** | 内核跑 `Efield` 全网格，`EvolveEPML` 之后再叠加修正（RESEARCH2 §2.2） | EB 靠近 PML 时共形更新 + Yee 式 PML 修正双重作用，吸收可能失效；当前算例 EB 在内部未触发 |
| R7 | **J 项无对偶权重** | 内核直接写 `−c2·dt·μ0·j`，等价于 Yee 的源项 | `J ≠ 0`（加粒子）时源项与共形几何不一致，破坏电荷守恒；应改为对偶面积积分或显式 abort |
| R8 | **F 场（divE 清洁）未保留** | 新内核没有 `E += c²dt·grad(F)` 段（`EvolveE.cpp:203-228`） | ECT + 场清洁组合会静默丢修正项 |
| R9 | **checkpoint 未含三个新字段** | 只加到 `FieldTypes` 枚举，未加 `Fields.H` 的 ECT 字段列表 | 重启后 `A_y/Lx/Lz` 为初值 1.0 → 几何错误（不崩，但结果错） |
| R10 | **`Lx`/`Lz` 双重角色** | 既是 `Ey` 的分子权重，又是 `Ez`/`Ex` 的分母 | 稳定化时"只抬分母"会引入通量不一致（RESEARCH2 §3.3） |
| R11 | **2D 的 `IntVect` 只有 2 分量** | 代码按 2D 写（`Ey_nodal_flag` 等），3D 推广时需重新核对 | 3D 移植易错 |

### 4.3 独立于 C 侧的问题（不在本文范围，但会误导判断）

| # | 问题 | 证据 |
|---|---|---|
| X1 | **CW 激光源跨 grid 不稳定** | `max_grid_size < n_cell` 时 200 步内到 1e308；单 grid 稳定。机理未查 |
| X2 | 米氏算例 `max_step = 1e6` + `stop_time`，单 grid 时 160² 也只需约 950 步 | — |
| X3 | `inputs_mie_*` 已重新生成（`max_grid_size = n`），但 **TM（`by`）档尚未重跑** | 只有 TE 的收敛阶数据 |

---

## 5. 验证与复现

```bash
cd <repo>
cmake --build build_2d -j 8

SRC=benchmarks/cases/c02_0_pec
BIN=build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB

# L2a：规则域，常数几何与真实几何都应到机器精度
mkdir -p /tmp/c1/ect_real && cd /tmp/c1/ect_real
sed 's/^\(algo\.maxwell_solver *= *ect\)$/\1\nalgo.ect_conformal_E = 2/' \
    <repo>/$SRC/inputs_ect_2d_noeb > inputs
OMP_NUM_THREADS=8 <repo>/$BIN inputs > run.log 2>&1
python3 <repo>/$SRC/verify_l2a.py . /tmp/l2/yee_noeb

# L3 45° TE
mkdir -p /tmp/c1/l3_45te && cd /tmp/c1/l3_45te
sed 's/^\(algo\.maxwell_solver *= *ect\)$/\1\nalgo.ect_conformal_E = 2/' \
    <repo>/$SRC/inputs_2d_oblique_te_45deg > inputs
OMP_NUM_THREADS=8 <repo>/$BIN inputs > run.log 2>&1
python3 <repo>/$SRC/analyze_oblique.py . 45 te

# 米氏收敛阶（单 grid；TE 三档，两开关）
for tag in c0 c2; do for n in 40 80 160; do
  d=/tmp/mie_$tag/ey_$n; mkdir -p $d; cp $SRC/inputs_mie_ey_$n $d/inputs
  [ $tag = c2 ] && sed -i 's/^algo\.maxwell_solver = ect$/&\nalgo.ect_conformal_E = 2/' $d/inputs
  (cd $d && OMP_NUM_THREADS=8 ../../../$BIN inputs > run.log 2>&1)
done; done
python3 $SRC/analyze_mie_convergence.py /tmp/mie_c0 ey
python3 $SRC/analyze_mie_convergence.py /tmp/mie_c2 ey
```

---

## 6. 提交记录

| 提交 | 内容 |
|---|---|
| `4206d85` | C 侧实现（8 个源文件）+ `verify_l2a.py` 噪声地板 |
| `e6d352d` | `HANDOFF.md` 更正米氏归因 |
| `78b31c2` | 米氏根因（多 grid）+ 后处理两处修正 + 收敛阶基线 |
| `f1c6190` | 合并到 `master`，tag `c-side-2d-v1` |

工作区干净，`git status` 无输出。

---

## 7. 建议的下一步（按性价比）

1. **R1 + R2（半天）**：把 `Ey` 掩膜改成 `A_y <= 0`，加 `A_stab` 并做 `volfrac` 扫描。
   这是"让已写好的共形几何真正生效"的最低成本改动，也可能让 L3/米氏的数字明显变化。
2. **R3（1–1.5 天）**：level-set 几何，把近场阶数从 1.06 推到接近 2。
3. **X1（0.5 天）**：查 CW 激光源跨 grid 的机理，解除生产算例的分块限制。
4. **R7/R8/R9（半天）**：J 项 abort 保护、F 段保留、checkpoint 字段。
5. **R4（1.5 天）**：3D 推广；**R5（1 天）**：AMR per-level。
