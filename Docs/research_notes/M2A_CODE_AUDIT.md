# M2a 代码实现审查报告

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 审查对象：`<repo>` 中 ECT 求解器的 **2D XZ TE 模**实现（M2a-1b 的 B 侧）
> 对照基线：`~/WarpX_changed`（上游 clone，`a389dc54a`，development 2026-10-29）
> 审查日期：2026-10-31　|　构建：`build_2d`（2D XZ）/ `build`（3D），EB=ON、MPI、双精度
> 审查方式：源码逐函数阅读 + 与上游 `diff` + 复跑验证算例 + 新增分辨率/稳定性实验 + 文献/开源项目对照
>
> 结论速览：**B 侧（Faraday/环量）实现正确、可复现（L1/L2a/L2b 通过），但 E 侧共形 Ampere 尚未实现；
> TE 的反射精度与 L2b 判据仍有明确缺口；初始化时序修复有效，但覆盖面需补。**

---

## 0. 后续处理状态（2026-10-31 晚，审计之后）

> 本节由后续修订追加，记录审计问题清单（§6）的处置结果。原报告正文保持原样。

| 审计项 | 级别 | 处置 | 证据 |
|---|---|---|---|
| `verify_l2b.py` 判据/脚本不稳健（`dz` 硬编码、载波峰值判时） | P1 | **已修复** | 见 §0.1；64/128/256 全部 PASS |
| E 侧仍是掩膜 Yee（C 侧未实现） | P1 | **按用户指示暂缓** | 设计见 `C_SIDE_DESIGN.md` |
| TE 面 BCK-only | P2 | 保持现状 + 新增角度算例量化 | 60° 斜入射 `\|r+1\|=0.0434` |
| 2D `S_stab` 无文献、未做稳定性扫描 | P2 | 新增 45° 2000 步冒烟算例（非完整扫描） | `inputs_ect_2d_45deg` |
| 初始化修复只覆盖 external-E 路径 | P2 | **未改**（记录为已知限制） | 见 §3.6 原文 |
| 借用索引 / `S` 就地覆盖 / `ly` 保守 | P3 | **未改**（当前不可达或属设计取舍） | 见 §3.1/§3.3/§3.7 |

### 0.1 `verify_l2b.py` 的修订（P1）

旧脚本的三个问题与对应改法：

1. **`dz` 硬编码默认 0.25** → 现在从运行的 `warpx_used_inputs`（或 `inputs`）解析
   `amr.n_cell` 与 `geometry.prob_lo/hi`，并打印 `dz = ... (parsed from inputs: yes/no)`；
   `zc` 也从 `my_constants.zc` 解析。
2. **用绝对到达时间 `t_ref` 判定** → 改为**往返时间** `t_ref − t_inc`。探测的是载波峰值，
   其相位随 `dz` 漂移，但 `t_inc` 与 `t_ref` 同向漂移，相减后抵消（审计 §4.2 的 128²/256² 现象）。
3. **峰值比 `r` 当门禁** → 降级为信息性输出（它是载波峰值比，不是反射系数）。

修订后（解析往返 `2(z_probe−z_wall)/c = 10.007 ns`；壁内移一层 `−2dz/c`）：

| n | dz (m) | ECT 往返 (ns) | ECT 误差 (ns) | Yee 往返 (ns) | Yee 误差 (ns) | 判定 |
|---|---|---|---|---|---|---|
| 64 | 0.2500 | 10.378 | +0.371 | 8.255 | −1.752 | PASS |
| 128 | 0.1250 | 10.024 | +0.017 | 9.199 | −0.808 | PASS |
| 256 | 0.0625 | 10.083 | +0.076 | 9.670 | −0.336 | PASS |

→ 三条判据（往返误差 ≤0.5 ns、ECT 比 Yee 更接近解析、`dz` 自动解析正确）在三个分辨率上一致成立；
旧脚本在 128²/256² 上会因绝对时间偏早 0.57/0.69 ns 而误判 FAIL。

### 0.2 本次新增的算例（回应 §6 的 P2/P3 与 §7 建议 3/5）

| 文件 | 内容 | 判定 | 结果 |
|---|---|---|---|
| `run_l2b_scan.sh` | L2b 分辨率扫描 64/128/256 | 往返时间判据 | 三档全 PASS |
| `inputs_ect_2d_45deg` + `check_stability.py` | 45° 斜面 TE，2000 步 | 有限性 / 能量漂移 / 无增长 | PASS（`Emax/E0=1.000189`，末段漂移 `2.1e-4`） |
| `inputs_3d_te_pulse` + `analyze_3d_te.py` | 3D TE 脉冲（32³、64³） | 波速 / 幅值 / `Bx/Ey` | PASS；64³ 误差减半（`dt` 0.393→0.200 ns） |
| `inputs_2d_l3_60deg` + `analyze_l3.py` | 60° 斜入射反射 | 沿束平均 `\|r+1\| ≤ 0.05` | PASS（`r = −1.0434`） |
| `make_l3_angle.py` | 由 `inputs_2d_l3` 生成任意角度算例 | — | — |

> 15° 掠入射**不适用**峰值法：反射束与入射束几乎同向，两脉冲在探针处的时间间隔
> `(s − z_P)/c` 只有 1.3–3.6 ns，小于脉宽（w=1.5 → ~3 ns），波峰重叠。需改窄脉冲或扩大域，
> 本轮未做，故未保留该算例。

---

## 0. 结论摘要

| 维度 | 结论 | 证据 |
|---|---|---|
| 代码规模 | 相对上游改 **6 个文件**，`+196 / −48` 行 | `diff` 统计（§2） |
| 几何（2D per-unit-y） | **正确**：规则网格精确退化（`Sx={0,dz}`、`Sz={0,dx}`、`ly={0,1}`） | L1；`CASES_2026-10-31.md` C1–C6 |
| 环量 Rho | **正确**：XZ 下 `Rhox`/`Rhoz` 是 3D 公式的合法约化；规则域逐位等于 Yee | 本次复跑 L2a：201 步 `max\|d\|=0` |
| B 更新 | **正确**：`idim<3` 使 `Bz` 参与；三分支（extended/available/intruded）沿用上游 | `EvolveB.cpp:264` 起 |
| E 更新 | **未共形**：仍是掩膜 Yee（`EvolveE.cpp:92-94`）；C 侧缺失 | `grep EvolveECartesianECT` 无命中 |
| 初始化时序 | **已修复**：`EvolveECTRho` 移到 `AddExternalFields` 末尾（`WarpXInitData.cpp:987-1000`） | L2a/L2b 由 FAIL → PASS |
| TE 反射精度 | **受限**：64² 时 `r=−0.9188`（TM 为 `−0.999361`）；随分辨率收敛（128² `−0.9632`、256² `−0.9921`） | 本次新增分辨率扫描（§4.2） |
| 判据稳健性 | **有缺陷**：`verify_l2b.py` 的 `dz` 默认硬编码 0.25、用载波峰值判时；高分辨率下绝对时间判据会误判 | 脚本源码 + 128/256 结果 |
| 稳定化 | **BCK-only**：TE 面不做 ECT 借用；2D `S_stab` 为 3 项约化（非文献原文） | `EmbeddedBoundaryInit.cpp:451-465`；`WarpXFaceExtensions.cpp:96-133` |
| 稳定性（初步） | 45° 斜面 TE 跑 2000 步，探针值有限、无发散 | 本次实验 `/tmp/stab45` |
| 未覆盖 | L3 斜入射、C 侧、2D/3D 一致性、AMR、PML+TE、重启/Python 路径 | §6 |

**一句话**：这是一次**正确但只完成一半**的实现——B 侧（Faraday）已经与上游 ECT 框架对齐并通过逐位退化验证；
把它当作"2D TE 可用"之前，必须先补 C 侧（E 场对偶共形 Ampere）并修 L2b 判据，否则 8% 量级的反射误差无法归因。

---

## 1. 审查方法与范围

1. **代码**：对 6 个改动文件做 `diff -u`，逐函数阅读（`EmbeddedBoundaryInit.cpp`、`WarpXFaceExtensions.cpp`、
   `EvolveB.cpp`、`EvolveECTRho.cpp`、`WarpXFaceInfoBox.H`、`WarpXInitData.cpp`），并与上游 `EvolveBCartesianECT`
   三分支、`EvolveECTRho` 3D 公式、`EvolveE` 掩膜路径逐行对照。
2. **数值**：
   - 复跑 `benchmarks/cases/c02_0_pec/run_l2.sh`（L2a/L2b，四算例）；
   - **新增**：把 L2b 的对齐 PEC 算例在 128²、256² 上重跑（ECT 与 Yee 各一），观察 `r` 与往返时间随分辨率的变化；
   - **新增**：把平面转到 45°、跑 2000 步的 TE 稳定性测试；
   - **新增**：3D TE（Ey+Bx）脉冲，确认 3D ECT 在时序修复后仍正常。
3. **文献/开源**：WarpX 上游、AMReX、Meep、Gkeyll 源码；Dey–Mittra、Xiao–Liu、Benkler、Jurgens–Taflove、
   Pember、Colella、Berger–Giuliani、Farjadpour、Kottke 等论文；Taflove & Hagness、LeVeque 教科书。
4. **未做**：AMR 多层、PML+TE、2D/3D 逐点一致、长时间（>10⁴ 步）全局稳定性、C 侧（代码尚不存在）。

---

## 2. 代码改动清单（6 文件，+196/−48）

| # | 文件 | 函数 / 行号（fork） | 改动 | 审查评价 |
|---|---|---|---|---|
| 1 | `EmbeddedBoundaryInit.cpp` | `MarkUpdateECellsECT`（:233，XZ 分支） | Ey 掩膜改由 `ly`（节点开放度）判定 | ✅ 与 `ly` 的定义自洽 |
| 2 | 同上 | `MarkUpdateBCellsECT`（:299，XZ 分支） | Bx/Bz 掩膜改用 `face_areas[0]/[2]` | ✅ 与 Bx/Bz 的对偶面一致 |
| 3 | 同上 | `MarkExtensionCells`（:370；BCK 块 :451-465） | 去掉 XZ 对 `idim=0/2` 的跳过；2D `S_stab`；TE 面直接 BCK（`S→S_stab`，`flag=bck_stabilized`） | ✅ 结构正确；⚠️ 精度取舍见 §3.7 |
| 4 | 同上 | `ComputeEdgeLengths`（:498；`ly` 块 :564-590） | XZ 跳过 `idim=1` 的置零；由 `lx/lz` 四边判据填 `ly` | ✅；⚠️ 判据偏保守（§3.1） |
| 5 | 同上 | `ComputeFaceAreas`（:592；XZ 分支 :609-659） | `face_areas[0]/[2]` 由 `getAreaFrac()` 填（regular→1、covered→0、cut→fraction） | ✅；规则/覆盖 box 单独处理正确 |
| 6 | 同上 | `ScaleAreas`（:714；XZ :726） | 去掉 XZ 对 0/2 的跳过；`full_area=cell_size[(idim+1)%3]·cell_size[(idim+2)%3]` | ✅ 2D 下即 `dz`/`dx` |
| 7 | `WarpXFaceExtensions.cpp` | `ComputeSStab`（:84；XZ :96-133） | 新增 `dim=0/2` 的 2D 分支（3 项 `max`） | ✅ 形式上是 3D 的合理约化；⚠️ 无文献原文（§3.7） |
| 8 | 同上 | `CountExtFaces`（:32；XZ :45、:67） | XZ 统计 0..2；`ReduceIntSum(...,3)` | ✅ |
| 9 | 同上 | `shrink_borrowing`（:240；:261） | 循环扩到 3 | ✅ |
| 10 | `WarpXFaceInfoBox.H` | :56 | `int vecs_size = 0;` | ✅（消除未初始化） |
| 11 | `EvolveB.cpp` | `EvolveBCartesianECT`（:228；循环 :264） | `idim < AMREX_SPACEDIM` → `idim < 3` | ✅ 必要且充分 |
| 12 | `EvolveECTRho.cpp` | `EvolveRhoCartesianECT`（:68；XZ :120、:147） | XZ 下计算 `Rhox`/`Rhoz` | ✅ 见 §3.2 的公式核对 |
| 13 | `WarpXInitData.cpp` | `AddExternalFields`（:963；修复 :987-1000） | 在施加 E 边界后重算 `ECTRhofield` | ✅ 修复有效；⚠️ 覆盖面 §3.6 |
| 14 | 同上 | `InitLevelData`（:1336 起） | 删除早于 E 注入的 `EvolveECTRho`；保留 `Efield_cp` 分支（:1356） | ✅ |

> 上游既有、**不是**本次改动（避免误判改动面）：`MarkUpdateECellsECT` 的声明与调用（`EmbeddedBoundaryInit.H:75`、
> `WarpXInitData.cpp` 旧调用点）、`EvolveBCartesianECT` 的调度（`EvolveB.cpp:117`）、`init_borrowing` 的 0..2 循环。

---

## 3. 逐项审查

### 3.1 2D per-unit-y 几何

- **`edge_lengths[1]`（:564-590）**：把 3D 的 y 向棱退化为节点开放度（0/1）；判据是"节点周围四条边任一 fully covered → 0"。
  - ✅ 与 `MarkUpdateECellsECT` 的 Ey 判据一致；`CellSize[1]=1` 使 `ly` 无需缩放。
  - ⚠️ 该判据是**保守**的：节点只要挨着一条覆盖边就被判为 closed，可能把本可更新的 Ey 冻结，影响贴面精度。
    建议与"仅当节点本身在 PEC 内才 closed"的变体做对比（§7）。
- **`face_areas[0]/[2]`（:609-659）**：用 AMReX `getAreaFrac()`（2D 下是边的开放分数）填，regular/covered box 单独处理。
  - ✅ 与 2D 的 Bx（z 边）/Bz（x 边）对偶面一致；L1 实测 `Sx={0,dz}`、`Sz={0,dx}` 精确成立。
  - ✅ `ScaleAreas`（:726）用 `cell_size[1]=1` 正确得到 `dz`/`dx`。
- **`MarkExtensionCells`（:451-465）**：TE 面不进入扩展/借用，而是直接 BCK：`S<S_stab` 时 `S→S_stab`、
  `area_mod→S_stab`、`flag=bck_stabilized`。
  - ✅ 与上游 BCK 分支的语义一致（`WarpXFaceExtensions.cpp:158-159` 的注释也说明 BCK 是兜底且精度低于 ECT）。
  - ⚠️ `S` 被就地覆盖：`face_areas` 不再是真实开放面积；`EvolveECTRho` 的 `Rho=circ/S` 因此用放大后的面积
    （这正是 BCK 的欠松弛），但会丢失真实面积，后续若启用 ECT 借用/诊断需保留原值（§3.7）。

### 3.2 Rho（Faraday 环量）

`EvolveECTRho.cpp` 的 2D 分支：

```
Rhox = (Ey(j)·ly(j) − Ey(j+1)·ly(j+1)) / Sx          （:120）
Rhoz = (Ey(i+1)·ly(i+1) − Ey(i)·ly(i)) / Sz          （:147）
```

- ✅ 与 3D 公式在"∂y=0、`Ez` 项相消、`k→j`"下逐项一致；符号被 `B -= dt·Rho` 吸收（L2a 逐位验证）。
- ✅ 单位：`Ey·ly` = V/m，`/S`（m）→ V/m²，`dt·Rho` → T，正确。
- ✅ `EvolveECTRho` 在 `EvolveE` 末尾重算（`WarpXPushFieldsEM.cpp:1074-1083`），B 侧消费的是"刚更新完的 E"的环量。
- ⚠️ 对 BCK 面，`S` 已被放大（§3.1），`Rho` 是欠松弛值；这是方法本身，不是 bug，但需在文档中写明。

### 3.3 B 更新（Faraday）

- ✅ `idim<3`（:264）使 `Bz` 参与；`Bfield[2]`、`face_areas[2]`、`ECTRhofield[2]`、`Venl[2]`、`flag_info_cell[2]`、
  `borrowing[2]` 在 2D 都已分配（`WarpX.cpp:2683-2750`）。
- ✅ 三分支逻辑（:285-395）与上游一致：extended 用 `Venl/S_mod` 记账，available/bck 用 `B -= dt·Rho`，
  intruded 用 `Venl/S`。
- ⚠️ **潜在缺陷**：`idim=0/2` 的借用邻居映射仍是 3D-y 索引（`kp = k + vec(1)`，:305/314/337/346）；
  当前 TE 面 `flag_ext_face=0`（BCK-only），该分支不可达；一旦启用 ECT 借用就会读错邻居，必须同步修（§7）。

### 3.4 E 更新（Ampere）

- ❌ 仍是掩膜 Yee：`EvolveE.cpp:92-94` 对 ECT 调用同一个 `EvolveECartesian`，只在 :171/:182/:193 用
  `eb_update_E` 跳过被覆盖点；**没有** `EvolveECartesianECT`，也没有 `dual_edge_length`/`dual_face_area`。
- 后果：cut cell 上 E 的更新不是共形的，Faraday（B）共形而 Ampere（E）阶梯，界面处 E/B 不一致；
  这与 Dey–Mittra/Xiao–Liu 的"E 与 H 成对共形"要求不符（§5）。
- 证据：TE 对齐 PEC 的 `r=−0.9188`（TM 为 `−0.999361`），且在 64² 上误差远大于 TM；随分辨率收敛
  （128² `−0.9632`、256² `−0.9921`，§4.2），说明当前误差有相当部分来自分辨率/掩膜，但 E 侧缺失
  注定无法把 TE 做到 TM 的 0.1% 量级。

### 3.5 掩膜与"贴面 cut cell"

- `MarkUpdateECellsECT`（:233）：Ex/Ez 用 `lx/lz==0`，Ey 用 `ly==0`；`MarkUpdateBCellsECT`（:299）：Bx/Bz 用
  `Sx/Sz==0`，By 用 `Sy==0`。
- **贴面 cut cell**：对齐 PEC 平面时，AMReX 把贴面单元标为 cut（`volfrac=1` 但 `isRegular()=false`）；
  Yee 的 `MarkUpdateCellsStairCase` 把它当封闭（多冻结一层），ECT 的几何（`getAreaFrac`/`edge_cent`）把它当开放。
  这是 L2b 差异的全部来源（报告 §14.3），也是"ECT 壁面在 z=0、Yee 在 z=−dz"的原因。
- ✅ 该行为与文献中"staircase vs cut-cell"的差异一致（Pember 1995；Colella 2006）；**不是 fork 的 bug**，
  但任何与 Yee 的逐位比较都必须先排除这一层约定差异。

### 3.6 初始化时序修复

- **缺陷**：旧代码在 `InitLevelData` 中、初始 E 注入 `Efield_fp` **之前**调用 `EvolveECTRho` → `Rho≡0` →
  第一个半步 `B -= dt·0` 空操作 → 整体滞后半步（报告 §7）。
- **修复**：删除旧调用；在 `AddExternalFields` 末尾（`ApplyEfieldBoundary` 之后）重算（`WarpXInitData.cpp:987-1000`）。
  - ✅ 同时覆盖初始化与演化中改场（`WarpXEvolve.cpp:330`）；对 constant 与 parser 两种外部 E 都生效；
    保留 `Efield_cp` 分支（:1356，其填场在同函数内先行完成）。
  - ✅ 重启：`ECTRhofield`/`Venl` 在 checkpoint 字段列表中（`Fields.H:155-156`），重启时随场恢复，通常无需重算。
  - ⚠️ **覆盖面**：修复只在 `E_ext_grid_type != default_zero` 分支内执行；若 E 由 Python 回调/其他初始化路径
    在 `AddExternalFields` 之后修改，`Rho` 仍可能陈旧。建议把"E 改动后重算 Rho"做成一个显式的
    `SyncECTRho(lev)`，在 `AddExternalFields`、重启、Python 回调之后统一调用（§7）。
  - ⚠️ 多级：ECT 几何只在 `maxLevel()` 初始化（`WarpXInitData.cpp:1607-1633`）；`AddExternalFields` 现在对
    每层都尝试重算，低层 `S=0` → `EvolveECTRho` 直接返回（安全但低层 Rho 仍为 0）。多层 ECT 本来就未支持，
    这里只需不让它崩溃/误导。

### 3.7 稳定化：BCK-only 与 2D `S_stab`

- **BCK-only**（:451-465）：TE 面（Bx/Bz）不做 ECT 借用，只在 `S<S_stab` 时把面积抬到 `S_stab` 并标
  `bck_stabilized`；B 更新走 `B -= dt·Rho`（不带 `S_mod`）。
  - ✅ 与上游 BCK 分支一致；`WarpXFaceExtensions.cpp:158-159` 的注释明确"BCK 精度低于常规 ECT，但优于阶梯"。
  - ⚠️ 因此 TE 的精度上限低于"完整 ECT"；若 L3 需要 <1% 反射误差，建议把 TE 面也做成 ECT 借用（同时修
    借用索引与 `Venl` 记账），或至少做 BCK vs ECT 的对比算例。
- **2D `S_stab`**（`WarpXFaceExtensions.cpp:96-133`）：`dim=0: 0.5·max{ly(j)·dz, ly(j+1)·dz, lz·dy}`；
  `dim=2` 对称。
  - ✅ 形式上与 3D 的 4 项公式在"y 退化"下一致（两个 y 边 + 一个 z 边）。
  - ⚠️ **文献中没有这个 2D 约化**（Xiao–Liu 2008 只给 3D 判据）；其稳定性/精度只能靠数值验证。
    本次 45° 测试只证明"2000 步不发散"，还不是稳定性边界（需角度/分辨率/体积分数扫描，§7）。

### 3.8 边界、PML、AMR

- **PML**：`EvolveEPML`/`EvolveBPML` 对 ECT 走 Yee 分支（`EvolveEPML.cpp:83`、`EvolveBPML.cpp:68`）；
  L2 算例用 PEC 边界，**TE+PML 未测**；TM 基线（PML 四边）正常（`r=−0.999361`）。
- **AMR**：几何只在 `maxLevel()` 初始化/扩展（`WarpXInitData.cpp:1607-1633`）；`face_areas` 的 per-unit-y
  约化在 coarse/fine 界面、`area_mod`/`Venl` 上未验证。M2a 目前固定单层，与计划一致，但需在 M4 前补。
- **2D/3D 一致性**：决策要求"2D TE 与 3D 逐点一致"，尚未搭建对照（3D 需 1 层 y + 周期 + 同 `dx/dz/dt`）。

---

## 4. 数值验证（本次复跑与新增实验）

### 4.1 L2a：规则域逐位退化 —— PASS

`run_l2.sh` 复跑（ECT vs Yee，EB 面移出域外，201 步）：

```
--- fp1 / fp2 / FieldEnergy (各 201 采样) ---
Ey  peak=6.235032e+05  max|d|=0.000000e+00
Bx  peak=1.944326e-03  max|d|=0.000000e+00
total peak=4.060209e+01 max|d|=0.000000e+00
RESULT: PASS
```

→ 在无 cut cell 时，per-unit-y 约化在**算术层面**与 Yee 完全等价；此前的 0.55% 偏差确系初始化时序缺陷。

### 4.2 L2b：对齐 PEC 的反射 —— 64² PASS，但判据不稳健

| n | dz (m) | ECT `r` | ECT 往返时间 (ns) | Yee `r` | Yee 往返时间 (ns) | 解析往返 (ns) |
|---|---|---|---|---|---|---|
| 64 | 0.25 | −0.9188 | 10.378 | −0.8505 | 8.255 | 10.007（壁在 0）/ 8.340（壁在 −dz） |
| 128 | 0.125 | −0.9632 | 10.024 | −0.9668 | 9.198 | 10.007 / 9.173 |
| 256 | 0.0625 | −0.9921 | 10.083 | −0.9923 | 9.670 | 10.007 / 9.590 |

- **往返时间**（`t_ref−t_inc`，消掉载波相位偏差）随分辨率收敛到解析值；ECT 在 ~0.2–0.8%（128/256），
  Yee 落在"壁内移一层"的解析值附近（~0.3–1.0%）→ L2b 的物理结论成立。
- **绝对时间 `t_ref`**：64² 为 21.700 ns（解析 21.682），但 128²/256² 为 21.110/20.992（偏早 0.57/0.69 ns）。
  原因是 `verify_l2b.py` 用**载波峰值**判时，而载波相位随 `dz`/采样变化；`t_inc` 也同样偏早
  （11.322 → 11.086 → 10.909），所以绝对时间不是稳健量。
- **脚本缺陷**：`verify_l2b.py:82` 的 `dz` 默认硬编码 `0.25`；在 128²/256² 上不传 `dz` 会得到错误的解析值
  （我复跑时确认）。此外"`|t_ref−解析|≤0.5 ns`"在 128/256 上按绝对时间会 FAIL，而按往返时间则 PASS。
- **反射幅值** `|r|` 随分辨率从 0.92 收敛到 0.99（256²），但这是**峰值比**而非严格的反射系数；需要
  双探针分解或 NTFF 才能给出可信的 `r`（尤其要和 TM 的 0.07% 比）。

### 4.3 45° 斜面稳定性 —— 初步通过

把 L2b 的平面转到 45°（`eb_implicit_function = cos(π/4)x + sin(π/4)z`），2000 步：

```
fp1 Ey: max=7.254e+05  min=-7.765e+05  finite=True
```

→ 无发散；但还未做角度/分辨率/体积分数扫描，不能当作"BCK-only 稳定"的证明。

### 4.4 3D TE —— 正常

3D 64³、TE（Ey+Bx）脉冲，时序修复后：质心 −2.03 → +0.90（t=9.63 ns，`c·t=2.887`，预期 +0.86）；
Ey/Bx 幅值稳定。→ 3D ECT 未被 2D 改动影响，且时序修复对 3D（初场 E≠0）同样有效。

### 4.5 其它回归（引用报告 §14.4）

| 算例 | 结果 |
|---|---|
| 2D `embedded_boundary_rotated_cube`（上游 TM 回归） | 相对差 6.8e-10 |
| 3D `embedded_boundary_rotated_cube` | ≤1e-13 |
| 2D TM 基线 `inputs_2d_tm_flat`（n=384） | `r=−0.999361`（修复前 −0.999295） |

---

## 5. 与开源项目 / 文献对照

| 主题 | 文献/实现 | 本 fork 的对应 | 评价 |
|---|---|---|---|
| 共形 FDTD（PEC） | Jurgens & Taflove 1992（contour path）；Dey & Mittra 1997（local conformal，**E 与 H 成对**） | 只把 Faraday（B）共形化；Ampere（E）仍掩膜 Yee | **主要缺口**：与 DM 的"成对共形"不一致；TE 精度受限于此 |
| ECT（面积扩展/借用） | Xiao & Liu 2004（enlarged cells）、2008（3-D ECT + 稳定性判据）；WarpX 上游实现（三分支、`Venl`、`area_mod`） | B 侧沿用上游；TE 面只做 BCK（Benkler 2006），不做借用 | **可接受但精度较低**；上游注释也承认 BCK < ECT |
| BCK 稳定化 | Benkler, Chavannes & Kuster 2006（几何精度 + 稳定性判据） | `MarkExtensionCells` 内直接抬面积（:451-465） | 与上游 BCK 分支一致；精度/稳定性需扫描 |
| 2D `S_stab` 约化 | Xiao–Liu 2008 只给 3D 判据 | 3 项 `max` 的 2D 版（`:96-133`） | **无文献原文**；需数值验证（本次仅 45° 初步） |
| 阶梯 vs cut-cell | Pember et al. 1995；Colella et al. 2006 | 贴面 cut cell 上 Yee 多冻结一层（L2b） | 与文献一致；是"约定差异"而非 bug |
| EB 有限体积/SRD | Berger & Giuliani 2021/2022/2024 | 尚未接入（M2a-2） | 计划中 |
| 介质界面（subpixel） | Farjadpour et al. 2006；Kottke et al. 2008；Deinega & Valuev 2007 | 尚未接入（M2b） | 计划中 |
| 矩方程 | Gkeyll（`moments/zero`，本地 `~/gkeyll`） | `Source/MomentsSolver` 已移植（F1–F3） | 与 M2a 无关；M2a-2/M2c 会用到 |
| FDTD 基础 | Taflove & Hagness（教科书）；LeVeque（FV 教科书） | — | 用于稳定性/色散判断 |
| 参考实现 | WarpX 上游 `a389dc54a`；AMReX（EB API/SRD）；Meep（subpixel） | 本 fork 改动见 §2 | 上游是唯一"逐位对照"基线 |

**关键对照结论**：

1. Dey–Mittra 的原始格式**同时**修改 Ampere 与 Faraday；本 fork 只修改了 Faraday。文献上这会导致界面处
   的 E/B 不一致，TE 反射精度应显著低于 TM（实测 8% vs 0.07%，与之一致）。
2. ECT 的"借用"与 BCK 的"抬面积"是两种不同稳定化；前者精度更高。本 fork 对 TE 选了 BCK-only，
   属于"先稳定、后精度"的工程折中，但必须在 L3 前决定是否升级。
3. 2D per-unit-y 约化在规则域上与 Yee 逐位一致（L2a），这是本次实现最扎实的部分；cut-cell 部分的行为
   仍需 C 侧与 L3 才能闭环。

---

## 6. 问题清单（按严重度）

| 级别 | 问题 | 影响 | 证据 |
|---|---|---|---|
| **P1** | E 侧仍是掩膜 Yee（C 侧未实现） | TE 精度上限低；Faraday/Ampere 不一致；L3 无法闭环 | `EvolveE.cpp:92-94`；`grep EvolveECartesianECT` 无命中 |
| **P1** | L2b 判据/脚本不稳健（`dz` 硬编码默认、载波峰值判时） | 高分辨率会误判；PASS 有偶然性；无法作为回归门禁 | `verify_l2b.py:82`；§4.2 的 128/256 数据 |
| **P2** | TE 面 BCK-only（不借用） | 精度低于 ECT；当前 64² `r=−0.9188` | `EmbeddedBoundaryInit.cpp:451-465`；上游注释 |
| **P2** | 2D `S_stab` 无文献原文、未做稳定性扫描 | BCK 阈值可能过强/过弱 | `WarpXFaceExtensions.cpp:96-133`；仅 45° 一例 |
| **P2** | 初始化修复只覆盖 external-E 路径 | Python/其他改场路径可能留下陈旧 Rho | `WarpXInitData.cpp:987-1000`；`WarpXEvolve.cpp:330` |
| **P2** | PML+TE、AMR、2D/3D 一致性未测 | 生产使用前必须补 | §3.8 |
| **P3** | `EvolveB` 的 `idim=0/2` 借用索引仍是 3D-y（`k+vec(1)`） | 当前不可达；启用借用时必错 | `EvolveB.cpp:305/314/337/346` |
| **P3** | `S` 被就地覆盖（BCK） | 丢失真实面积，阻碍诊断/后续借用 | `EmbeddedBoundaryInit.cpp:459` |
| **P3** | `ly` 的节点判据偏保守 | 贴面 Ey 可能被过度冻结 | `EmbeddedBoundaryInit.cpp:564-590` |

---

## 7. 建议（按优先级）

1. **修 L2b 判据与脚本**（当天可做）：
   - `verify_l2b.py` 从 `warpx_used_inputs`/`inputs` 解析 `dz`，去掉硬编码默认；
   - 改用**往返时间**（`t_ref−t_inc` vs `2L/c`）或**包络峰**，而不是载波峰值；
   - 反射系数改用双探针分解（入射/反射分离），不要用峰值比。
2. **C 侧（E 场对偶共形 Ampere）**（M2a-1b 的剩余部分）：
   - 按 `M2A_SOURCE_REVIEW_FINDINGS.md` §5.2 的 2D 语义建 `dual_edge_length`/`dual_face_area`；
   - E 侧同样要有 `S_stab`/借用/BCK 兜底；与 B 侧的 `Venl` 记账配对；
   - 先做"退化回 Yee"单测，再接 L3。
3. **TE 精度决策**（C 侧之后）：
   - 对比 BCK-only 与 ECT 借用（同时修 `EvolveB` 借用索引）；
   - 做角度（15°/30°/45°/60°）× 分辨率 × 体积分数扫描，给出 `r` 的误差带；
   - 目标：TE 反射误差 ≤1%（向 TM 的 0.07% 看齐）。
4. **初始化覆盖**：
   - 把"E 改动后重算 Rho"抽成 `SyncECTRho(lev)`，在 `AddExternalFields`、重启、Python 回调后统一调用；
   - 加重启算例（checkpoint→restart，比较场与 Rho）。
5. **回归与测试**：
   - 几何单测（平网格=1、半覆盖=0.5、错一格不串值）；
   - PML+TE、2D/3D 逐点一致（1 层 y + 周期）、AMR 单层→两层；
   - 长时间（≥10⁴ 步）稳定性/能量扫描。
6. **文档**：把"贴面 cut cell 的 ECT/Yee 约定差异"写进用户文档，避免把 L2b 的 8% 当 bug。

---

## 8. 参考文献（DOI 已核对）

1. K. S. Jurgens, A. Taflove, *Finite-difference time-domain modeling of curved surfaces*,
   IEEE Trans. Antennas Propag. 40(4):357–366 (1992). DOI: 10.1109/8.138836
2. S. Dey, R. Mittra, *A locally conformal finite-difference time-domain (FDTD) algorithm for modeling
   three-dimensional perfectly conducting objects*, IEEE Microw. Guided Wave Lett. 7(9):273–275 (1997).
   DOI: 10.1109/75.622536
3. T. Xiao, Q. H. Liu, *Enlarged cells for the conformal FDTD method to avoid the time step reduction*,
   IEEE Microw. Wirel. Components Lett. 14(12):551–553 (2004). DOI: 10.1109/LMWC.2004.837384
4. T. Xiao, Q. H. Liu, *A 3-D enlarged cell technique (ECT) for the conformal FDTD method*,
   IEEE Trans. Antennas Propag. 56(3):765–773 (2008). DOI: 10.1109/TAP.2008.916876
5. S. Benkler, N. Chavannes, N. Kuster, *A new 3-D conformal PEC FDTD scheme with user-defined geometric
   precision and derived stability criterion*, IEEE Trans. Antennas Propag. 54(6):1843–1849 (2006).
   DOI: 10.1109/TAP.2006.875909
6. R. B. Pember, J. B. Bell, P. Colella, W. Y. Crutchfield, M. L. Welcome, *An adaptive Cartesian grid
   method for unsteady compressible flow in irregular regions*, J. Comput. Phys. 120(2):278–304 (1995).
   DOI: 10.1006/jcph.1995.1165
7. P. Colella, D. T. Graves, B. J. Keen, D. Modiano, *A Cartesian grid embedded boundary method for
   hyperbolic conservation laws*, J. Comput. Phys. 211(1):347–366 (2006). DOI: 10.1016/j.jcp.2005.05.026
8. M. J. Berger, A. Giuliani, *A state redistribution algorithm for finite volume schemes on cut cell
   meshes*, J. Comput. Phys. 428:109820 (2021). DOI: 10.1016/j.jcp.2020.109820
9. A. Giuliani et al., *A weighted state redistribution algorithm for embedded boundary grids*,
   J. Comput. Phys. 464:111305 (2022). DOI: 10.1016/j.jcp.2022.111305
10. M. J. Berger, A. Giuliani, *A new provably stable weighted state redistribution algorithm*,
    SIAM J. Sci. Comput. 46(5):A2848–A2873 (2024). DOI: 10.1137/23M1597484
11. A. Farjadpour et al., *Improving accuracy by subpixel smoothing in the finite-difference time domain*,
    Opt. Lett. 31(20):2972–2974 (2006). DOI: 10.1364/OL.31.002972
12. A. Kottke, A. Farjadpour, S. G. Johnson, *Perturbation theory for anisotropic dielectric interfaces...*,
    Phys. Rev. E 77:036611 (2008). DOI: 10.1103/PhysRevE.77.036611
13. W. Zhang et al., *AMReX: a framework for block-structured adaptive mesh refinement*,
    J. Open Source Softw. 4(37):1370 (2019). DOI: 10.21105/joss.01370
14. A. Taflove, S. C. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain Method*,
    3rd ed., Artech House (2005).
15. R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge Univ. Press (2002).

**开源项目**：WarpX 上游（`github.com/ECP-WarpX/WarpX`，本机 `~/WarpX_changed`）；
AMReX（`~/amrex-development`）；Meep（`<local-ref>/meep`）；
Gkeyll（`~/gkeyll`）。

---

## 附录 A　复现命令与本次结果

```bash
# L2a/L2b（四算例，约几分钟）
cd <repo>/benchmarks/cases/c02_0_pec && ./run_l2.sh

# 本次新增：分辨率扫描（128²/256²）
#   用同一 inputs_ect_2d_degen / inputs_yee_2d_degen，改 amr.n_cell 与 max_step
#   结果见 §4.2；数据在 /tmp/scan/{ect,yee}{128,256}

# 本次新增：45° 稳定性（2000 步）
#   eb_implicit_function = cos(pi/4)x + sin(pi/4)z；结果 /tmp/stab45

# 本次新增：3D TE 复核
#   /tmp/ect_test/r3d_te（64³）；结果见 §4.4
```

关键数字：L2a `max|d|=0`（201 步）；L2b 64² `r=−0.9188`（ECT）/`−0.8505`（Yee）；128² `−0.9632`/`−0.9668`；
256² `−0.9921`/`−0.9923`；45° 2000 步 `|Ey|max=7.8e5`、有限；3D TE 质心 −2.03→+0.90。

## 附录 B　行号索引（fork 当前工作树）

| 主题 | 位置 |
|---|---|
| E 掩膜（XZ） | `EmbeddedBoundaryInit.cpp:233-291` |
| B 掩膜（XZ） | `:299-368` |
| 扩展分类 + BCK（XZ） | `:370-471`（BCK 块 :451-465） |
| `edge_lengths[1]`（XZ） | `:498-590`（ly 块 :564-590） |
| `face_areas[0]/[2]`（XZ） | `:592-659`（getAreaFrac :609/614） |
| `ScaleAreas`（XZ） | `:714-732` |
| `ComputeSStab`（XZ） | `WarpXFaceExtensions.cpp:96-133` |
| `CountExtFaces`（XZ） | `:40-70` |
| `shrink_borrowing` | `:258-261` |
| B 三分支 / `idim<3` | `EvolveB.cpp:228`、`:264`、`:285-395` |
| 借用索引（潜在） | `EvolveB.cpp:305/314/337/346` |
| Rho（XZ） | `EvolveECTRho.cpp:120`、`:147` |
| E 更新（共用 Yee） | `EvolveE.cpp:92-94`、`:171/182/193` |
| 时序修复 | `WarpXInitData.cpp:987-1000`；保留的 cp 分支 `:1356` |
| checkpoint 字段 | `Fields.H:155-156`（`ECTRhofield`、`Venl`） |

## 附录 C　上游 ECT 与 fork 的对应关系

| 上游行为 | fork 的 2D 扩展 |
|---|---|
| `EvolveBCartesianECT` 循环 `idim<AMREX_SPACEDIM` | 改 `idim<3`（Bz 参与） |
| `EvolveRhoCartesianECT` 在 XZ 只算 `Rhoy` | 新增 `Rhox`/`Rhoz`（Ey 项） |
| XZ 下 `face_areas[0]/[2]=0`、`edge_lengths[1]=0` | 由 `getAreaFrac`/节点开放度填充 |
| XZ 只扩展 `idim=1` | TE 面改为 BCK-only（不扩展） |
| `ComputeSStab` 只有 3D 分支 | 新增 2D `dim=0/2` 分支 |
| 初始化只在 `InitLevelData` 里算一次 Rho（早于 E 注入） | 移到 `AddExternalFields` 末尾重算 |

> 本报告只审查 M2a 的 B 侧实现；C 侧（E 场对偶共形 Ampere）尚未写代码，不在审查范围内。
</content>
