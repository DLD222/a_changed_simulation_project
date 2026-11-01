# M4（AMR）实现调研与建议

> 面向 `WarpX_fluid` 的分层自适应网格：文献调研 → 代码落点 → 实施建议
>
> - 范围：M4 = 每层 EB/材料分配与初始化、coarse–fine 时空插值、reflux、子循环策略、regrid 重建
> - 状态：**调研完成，未实现**。本文不含任何仿真结果；所有"文献结论"均给出 DOI 或版次。
> - 依据代码：`<repo>`（WarpX 26.10 分支 + `Source/MomentsSolver`），
>   AMReX 源码：`~/amrex-development`
> - 写作日期：2026-10-31

---

## 0. 结论摘要（先看这 8 条建议）

| # | 建议 | 关键依据 |
|---|------|---------|
| R1 | **用 AMReX `FillPatcher` 做 coarse–fine 的时空插值**，不要自己写 | AMReX 提供 `AMReX_FillPatcher.H`；WarpX 全仓库对它 **0 引用**（`grep -rn FillPatcher Source/` 无输出），说明这是"未使用"而非"被删除" |
| R2 | **解除 `lev==maxLevel()` 的 EB/材料限制**：EB 几何、`edge_lengths/face_areas/area_mod`、材料系数、`eb_update` 都按层分配 | 当前 `WarpX.cpp:2631`（材料禁 AMR）、`WarpXInitData.cpp:1607-1632`（ECT 只在 maxLevel）；EB+AMR 的按层几何在文献中有成熟做法 [11] |
| R3 | **第一版强制 `do_subcycling=0`**（两层同 dt），把子循环留到最后 | 子循环下的界面守恒需要 flux correction [2][3]；先去掉一个自由度 |
| R4 | **reflux 分两类**：场通量用 Berger–Colella 通量校正 [2]；粒子/流体电流用"电荷–电流校正" [17][19] | WarpX 已有 `SyncCurrentAndRho` 的框架可复用 [18] |
| R5 | **材料系数的粗化必须是守恒的**（体积分数加权），界面处用 subpixel 平均 | [14][15]；WarpX 现有 `MacroscopicProperties` 是 cell-centered 常数/parser，未做 per-level |
| R6 | **椭圆/隐式求解直接用 AMReX MLMG**（多层几何多重网格），不自研 | [6] 给出 cell-centered AMR 投影法的多层椭圆求解；AMReX 提供 MLMG [4][5] |
| R7 | **regrid 钩子里重建 EB/材料/几何缓存**（`area_mod`、borrowing、`flag_info_face`） | [11] 的 per-level geometry；WarpX 已有 `build_eb_data_per_level` 开关（`WarpXInitEB.cpp:88`） |
| R8 | **验收用 C10（核壳球 + AMR）**，并新增三条 AMR 专属判据：coarse–fine 无伪反射、守恒、收敛阶 | [2] 的守恒测试范式；C10 见配套文档 |

---

## 1. 现状：本机代码在 AMR 上的真实边界

### 1.1 已经有的（可直接用）

| 能力 | 位置 | 说明 |
|---|---|---|
| 块结构 AMR 骨架 | AMReX `AmrCore` | `WarpX` 继承 `amrex::AmrCore` |
| 场边界填充 | `WarpX::FillBoundaryE/B` | 只做同层/邻层 `FillBoundary`，**不做时间插值** |
| 粗化 | `AverageDown*` | 场量守恒粗化 |
| 电流/电荷同步 | `SyncCurrentAndRho` | 粒子电流在 coarse–fine 界面的校正 |
| EB 每层几何开关 | `warpx.build_eb_data_per_level`（`WarpXInitEB.cpp:88`） | 但场求解器侧只在 maxLevel 用 |
| 多层椭圆求解 | AMReX `MLMG` | WarpX 的静电/隐式已用 |

### 1.2 明确缺失/受限的（M4 要补的）

| 限制 | 位置 | 后果 |
|---|---|---|
| 宏观材料禁止 AMR | `WarpX.cpp:2631-2632` | 材料只能在单层 |
| ECT 度量只在 `maxLevel` | `WarpXInitData.cpp:1607-1632` | 多层 EB cut-cell 无几何 |
| `eb_update_E/B` 只在 `maxLevel` 构造 | `WarpX.cpp`（分配段） | 低层场更新会空指针 |
| 没有 coarse–fine 时间插值 | 全仓库无 `FillPatcher` 引用 | 两层不同 dt 时场在界面上不闭合 |
| MomentsSolver 无 AMR | `Source/MomentsSolver/FV/MomentFVDriver.H` | 1D、单层、AMReX-independent |
| 混合/Darwin 单层 | `WarpXPushFieldsHybridPIC.cpp`（`finest_level==0`） | 与 AMR 互斥 |

> 注意：`FillPatcher` 零命中**不等于**"基础设施被移除"——AMReX 里有
> `Src/AmrCore/AMReX_FillPatcher.H`（类 `FillPatcher`，构造参数为 fine/coarse 的
> BoxArray、DistributionMapping、Geometry），只是 WarpX 没用它。上一份 `FINDINGS.md`
> 把"零命中"解释成"被裁掉"是**缺少基线对照**的误判。

---

## 2. 文献综述（按主题：文献说了什么 → 对本项目意味着什么）

### 2.1 块结构 AMR 的骨架

- **Berger & Oliger (1984)** 建立了块结构自适应网格的基本框架：按误差估计加密矩形块，
  块之间用插值耦合 [1]。
- **Berger & Colella (1989)** 给出守恒型 AMR 的完整算法（含 reflux），是后来几乎所有
  AMR 流体/PIC 代码的模板 [2]。
- **AMReX**（JOSS 2019 [4]；IJHPCA 2021 [5]）把上述算法工程化：`BoxArray`、`DistributionMapping`、
  `AmrCore`、`FillBoundary`、`AverageDown`、`MLMG`、`FillPatcher`、EB 工具。

**对本项目的含义**：M4 不应从头写 AMR；应把 `MomentsSolver` 的状态放进 AMReX 的
`MultiFab`/`AmrCore` 体系，复用其网格管理与椭圆求解。AMR 的算法与工程综述见
Plewa–Linde–Weirs 的专著 [22]。

### 2.2 子循环（subcycling）与界面守恒

- **Berger (1987)** 证明：细层比粗层多走几步时，界面上的通量不平衡会破坏守恒；必须做
  **flux correction（reflux）** [3]。
- **Berger & Colella (1989)** 给出了 reflux 的标准实现：把细层界面通量对粗层做累加修正 [2]。

**含义**：`warpx.do_subcycling=1` 不能在"没有 reflux"的情况下使用。对流体矩方程，
reflux 要把守恒量（密度、动量、能量/压力张量）的通量都算进去；对粒子，见 §2.3。

### 2.3 PIC 的电荷/电流守恒与 AMR

- **Esirkepov (2001)** 给出任意形状因子的精确电荷守恒电流沉积 [19]；这是 ∇·E=ρ 的离散基础。
- **Vay et al. (2002)** 系统讨论了 PIC 的 mesh refinement：细层边界上的电流/电荷需要
  **校正**，否则会出现自力和电荷不守恒 [17]。
- **WarpX (2021)** 在生产规模上实现了 AMR + 多种求解器，报告了 coarse–fine 的处理方式 [18]。

**含义**：粒子路径的 reflux 可复用 WarpX 现成的 `SyncCurrentAndRho` 思路；
流体路径（矩方程）没有现成实现，需要按 §2.2 的通量校正自己写，
并且要与粒子路径的电荷/电流校正**统一**（否则混合模型里 ∇·E−ρ 会漂移）。

### 2.4 多层椭圆求解

- **Martin & Colella (2000)** 给出 cell-centered AMR 的投影法：多层 Poisson 用
  coarse–fine 边界条件 + 多层 V-cycle [6]。
- **AMReX MLMG** 是这套算法的工程实现 [4][5]。

**含义**：Poisson（静电、散度清理、隐式 Maxwell）直接用 MLMG；不要自研多层求解器。
需要注意：MLMG 要求每层的 `BoxArray`/`DistributionMapping` 与场一致，且 coarse–fine
边界条件类型（Dirichlet/Neumann）要与 EB 边界一致。

### 2.5 EB 与 AMR

- **Pember et al. (1995)** 给出 AMR + Cartesian embedded boundary 的守恒算法：
  每层独立构建几何，粗化/插值要带几何信息 [11]。
- **Osher & Sethian (1988)** 的 level set 是表示任意几何的标准工具 [12]。
- **Dey & Mittra (1997)** 的 conformal FDTD 给 PEC 边界的 cut-cell 离散 [13]。

**含义**：EB 几何必须**每层独立**构建，且 coarse–fine 界面上的 cut-cell 面元要一致；
`area_mod`/borrowing 的扩展在粗化后不能破坏守恒。这是 M4 最容易被低估的部分。

### 2.6 材料界面的亚网格平均

- **Farjadpour et al. (2006)** 的 subpixel smoothing：用体积/面积分数与法向投影得到
  有效的 1/ε，显著优于阶梯，接近二阶 [14]。
- **Deinega & Valuev (2007)** 把 subpixel 推广到有耗/色散介质 [15]。
- **Xiao & Liu (2004)** 的 enlarged cells：在 conformal FDTD 里扩大过小的面元以避免
  CFL 缩水（WarpX 的 `area_mod` 即此思想）[16]。

**含义**：材料系数的"粗化"不能简单取算术平均；应当用体积分数加权的
（法向调和、切向算术）平均，并在 coarse–fine 界面上保证
`ε_coarse = (1/V)∫ε_fine dV` 的离散等价。否则界面会在粗层"漂移"。

### 2.7 流体/MHD 的 AMR 经验

- **Zachary, Malagola & Colella (1994)**：AMR + 高阶 Godunov 的 MHD 实现 [7]。
- **Balsara (2001)**：divergence-free AMR MHD——**prolongation 必须保持 ∇·B=0** [8]。
- **RAMSES (2002)** [9]、**Athena++ (2020)** [10]：生产级 AMR MHD 的设计与陷阱。

**含义**：矩方程若与 Maxwell 耦合（10 矩 + 磁场），粗化/插值必须处理 ∇·B=0；
重构、Riemann 求解与守恒的教科书参考见 LeVeque [23]；
可选 constrained transport 或散度清理。RAMSES/Athena++ 的经验表明：
**AMR 的精度瓶颈通常在 prolongation 与 reflux，而不是重构格式**。

### 2.8 矩/动理学流体

- **Wang, Hakim, Bhattacharjee, Germaschewski (2015)**：5 矩/10 矩与 PIC 的系统对比 [21]，
  给出闭合的选择依据。
- **Juno et al. (2018)**：Gkeyll 的 DG 动理学求解器 [20]，展示高阶方法在动理学问题上的收益。
- **Birdsall & Langdon**：PIC 噪声与色散的经典教材 [24]。

**含义**：矩方程的 AMR 验证必须同时对照 (i) 解析色散、(ii) PIC 基线；
AMR 本身不应改变色散关系（在均匀区），否则说明 prolongation 引入了伪源项。

---

## 3. 实施建议（每条：做什么、为什么、代码落点、风险）

### R1　用 AMReX `FillPatcher` 做 coarse–fine 时空插值

**做什么**：对 E/B、流体矩、材料系数，用 `amrex::FillPatcher` 在粗层时间步之间做
空间 + 时间的插值，填充细层 ghost；替代目前"只 `FillBoundary`、不做时间插值"的做法。

**为什么**：两层不同 dt 时，细层在中间时刻需要粗层的**时间插值**；Berger–Colella 的
算法要求 coarse–fine 边界上的场在时间上闭合 [2]。AMReX 已把该逻辑封装在
`AMReX_FillPatcher.H`（fine/coarse 的 BA、DM、Geometry 构造）[4][5]。

**代码落点**：`Source/FieldSolver/WarpXPushFieldsEM.cpp` 的 `FillBoundary*` 调用点；
新增 `WarpX::FillPatchFields(lev)`；在 `Evolve` 的每层推进前后调用。

**风险**：WarpX 的场是 staggered（nodal/face/edge），`FillPatcher` 需要正确的
`IndexType`；EB 场还需要每层几何一致（见 R2）。

### R2　解除 `lev==maxLevel()` 限制：EB/材料/度量按层分配

**做什么**：
1. `WarpX.cpp` 的材料分配（`MacroscopicProperties::AllocateLevelMFs`）从 `lev==0` 放开到每层；
2. `WarpXInitData.cpp` 的 ECT 分支从 `if (lev == maxLevel())` 改成每层；
3. `eb_update_E/B`、`edge_lengths/face_areas/area_mod`、`flag_info_face/ext_face` 每层分配；
4. `WarpXFaceExtensions` 的 borrowing 每层计算。

**为什么**：文献中的 EB+AMR 实现都是每层独立几何 [11]；当前"只在最细层"的写法在
两层时会让低层 `EvolveE` 解引用空指针（代码事实）。

**代码落点**：`Source/WarpX.cpp`（`AllocLevelMFs`）、`Source/Initialization/WarpXInitData.cpp`
（`InitializeEBGridData`）、`Source/EmbeddedBoundary/WarpXFaceExtensions.cpp`。

**风险**：内存与初始化时间随层数增长；`build_eb_data_per_level=1` 时 `m_eb_index_space`
的索引要按层维护（`WarpX.H:1089` 已有辅助函数）。

### R3　第一版强制 `do_subcycling=0`

**做什么**：M4 的第一版只支持"所有层同 dt"，把子循环作为后续增量。

**为什么**：子循环把 reflux、时间插值、稳定性三件事耦合在一起；Berger (1987) 说明
没有 reflux 的子循环不守恒 [3]。先同 dt 可以把 M4 拆成"几何/插值"与"时间"两件事。

**代码落点**：`WarpX.cpp` 已有 `do_subcycling` 的兼容性断言（`WarpX.cpp:765-778`），
M4a/M4b 阶段把材料/EB 与 `do_subcycling=1` 的组合显式禁止。

### R4　reflux：通量校正 + 电荷/电流校正

**做什么**：
- 场/流体：在 coarse–fine 界面上累加细层通量、修正粗层守恒量（Berger–Colella 的标准 reflux）[2]；
- 粒子：沿用/扩展 `SyncCurrentAndRho`，保证细层边界上的电流/电荷与粗层一致 [17][19]；
- 混合模型：粒子 + 流体电流合并后再做一次校正，保证 ∇·E−ρ 不漂移 [18]。

**为什么**：这是守恒的硬要求，不是精度优化 [2][3]。

**代码落点**：`Source/Parallelization/`（现有 `SyncCurrentAndRho` 的位置）；
流体侧新增 `reflux_fluid(lev)`。

**风险**：流体矩的守恒量包括压力张量（10 矩），reflux 要对所有分量做；
EB 面元上的通量要带几何因子，否则会与 cut-cell 的面积不一致（[11] 的教训）。

### R5　材料系数：守恒粗化 + subpixel

**做什么**：
1. 每层独立计算 ε/μ/σ（cell-centered），由 level set 采样得到；
2. 粗化时用体积分数加权：`ε_c = Σ f_i ε_i`（算术）与 `1/ε_c = Σ f_i / ε_i`
   （法向调和）分别用于不同分量 [14]；
3. 界面处用 Meep 式法向投影：`1/ε_eff,i = n_i² ⟨1/ε⟩ + (1−n_i²)/⟨ε⟩` [14]；
4. 有耗时用 [15] 的推广。

**为什么**：阶梯界面的误差是一阶的；subpixel 平均在曲面/斜面上接近二阶 [14][15]。
本项目的 C2（斜平面 cut-cell）就是用来量化这一步的。

**代码落点**：`MacroscopicProperties::AllocateLevelMFs` + 新的"per-level 材料初始化"；
与 EB 的 `area_mod`/`volume fraction` 对接。

**风险**：`ε` 的算术/调和平均在 coarse–fine 界面上可能不满足
`ε_coarse = average(ε_fine)`；必须用**同一套离散定义**在两层上计算，
否则界面会在粗层移动。

### R6　椭圆/隐式：用 AMReX MLMG

**做什么**：Poisson（静电、散度清理）、隐式 Maxwell 的多层求解走 MLMG。

**为什么**：[6] 的多层投影法就是 MLMG 的算法基础；AMReX 已实现并验证 [4][5]。

**代码落点**：`Source/FieldSolver/` 中已有 `MLMG` 的使用（静电/隐式），
扩展到"每层 + EB"的边界条件。

**风险**：MLMG 的 coarse–fine 边界条件与 EB 边界条件要在同一层上一致；
否则会出现"界面处电势不连续"的假象。

### R7　regrid 钩子：重建 EB/材料/几何缓存

**做什么**：在 `WarpX::regrid` 之后，按层重建 EB 几何、材料、`area_mod`、borrowing、
`flag_info_face`，并重新初始化 `eb_update`。

**为什么**：AMR 的网格在演化中会变；几何缓存必须跟着变 [11]。

**代码落点**：`Source/Parallelization/WarpXRegrid.cpp`；
`WarpX::InitEB`/`InitializeEBGridData` 的按层调用。

**风险**：regrid 期间的粒子/流体状态需要 remap；粒子有 `Redistribute`，
流体矩的状态需要 `average_down`/`interpolate`（§2.7 的 prolongation 问题）。

### R8　验证：C10 + 三条 AMR 专属判据

**做什么**：
1. **coarse–fine 无伪反射**：平面波穿过粗/细层界面，反射 < 1e-3；
2. **守恒**：总电荷、总能量、总动量在 reflux 前后闭合到 1e-6（相对）；
3. **收敛**：固定物理问题，`amr.max_level=0/1/2` 的结果差随细层分辨率下降。

**为什么**：这是 Berger–Colella 一类代码的标准验收 [2]；也是 C10 的判据来源。

**代码落点**：`benchmarks/cases/c10_coated_sphere/` + 新增 `c11_amr_wave/`（可选）。

---

## 4. 分阶段实施（建议顺序）

| 阶段 | 内容 | 验收 | 依赖 |
|---|---|---|---|
| M4a | 每层 EB/材料/度量分配与初始化；**不改时间推进** | 单层结果与 M2 逐位一致；两层静态几何无空指针 | R2, R7 |
| M4b | 两层、`do_subcycling=0`、`FillPatcher` 时空插值 | 平面波穿过 coarse–fine 界面，伪反射 <1e-3；C3 平板结果不变 | R1, R3 |
| M4c | 材料界面跨 coarse–fine；守恒粗化 + subpixel | C10a vs C10b 峰位差 ≤0.2%；能量/电荷守恒 | R4, R5 |
| M4d | 流体矩方程 + 粒子/流体耦合的 AMR | 均匀区色散不变；Sod 在 AMR 上收敛；C8/C9 基线不变 | R4, R6, R8 |

> 建议把 M4a/M4b 作为"几何与插值"的独立里程碑，M4c 才引入材料跨层，
> M4d 再动矩方程——这样每一步都有可判定的回归。

---

## 5. 风险与开放问题

| 风险 | 说明 | 缓解 |
|---|---|---|
| EB 几何跨层不一致 | 粗层与细层的 cut-cell 面元在界面上不匹配，破坏守恒 [11] | 每层独立几何 + 界面通量校正；C10 检验 |
| 材料界面跨 coarse–fine | ε 的平均在两层不一致会导致界面"漂移" [14] | 同一离散定义；体积分数加权；C10 峰位判据 |
| 子循环 + 矩方程 | 稳定性与 reflux 耦合 [2][3] | 先 `do_subcycling=0`；M4d 后再评估 |
| ∇·B=0 的 prolongation | 10 矩 + Maxwell 时粗化/插值破坏散度 [8] | constrained transport 或散度清理 |
| 性能 | 每层 EB/材料 + MLMG 的开销 | 先正确后优化；AMReX 的 GPU 后端可用 [5] |

**开放问题（需实验回答，不应在文档里假设）**：
1. WarpX 现有 `SyncCurrentAndRho` 是否足以覆盖"流体电流"的 coarse–fine 校正？
2. `FillPatcher` 对 staggered 场（E/B）的直接支持程度（需要试编译验证）；
3. 材料 subpixel 在 coarse–fine 界面上的"守恒粗化"是否与细层独立计算一致（C10 可测）。

---

## 6. 参考文献（DOI 均经 Crossref 核对，2026-10-31）

**AMR 基础**
1. M. J. Berger, J. Oliger, *Adaptive mesh refinement for hyperbolic partial differential
   equations*, J. Comput. Phys. 53(3), 484–512 (1984).
   DOI: [10.1016/0021-9991(84)90073-1](https://doi.org/10.1016/0021-9991(84)90073-1)
2. M. J. Berger, P. Colella, *Local adaptive mesh refinement for shock hydrodynamics*,
   J. Comput. Phys. 82(1), 64–84 (1989).
   DOI: [10.1016/0021-9991(89)90035-1](https://doi.org/10.1016/0021-9991(89)90035-1)
3. M. J. Berger, *On conservation at grid interfaces*, SIAM J. Numer. Anal. 24(5), 967–984 (1987).
   DOI: [10.1137/0724063](https://doi.org/10.1137/0724063)

**框架与求解器**
4. W. Zhang *et al.*, *AMReX: a framework for block-structured adaptive mesh refinement*,
   J. Open Source Softw. 4(37), 1370 (2019).
   DOI: [10.21105/joss.01370](https://doi.org/10.21105/joss.01370)
5. W. Zhang *et al.*, *AMReX: Block-structured adaptive mesh refinement for multiphysics
   applications*, Int. J. High Perform. Comput. Appl. 35(6), 512–526 (2021).
   DOI: [10.1177/10943420211022811](https://doi.org/10.1177/10943420211022811)
6. D. F. Martin, P. Colella, *A cell-centered adaptive projection method for the
   incompressible Euler equations*, J. Comput. Phys. 163(2), 271–312 (2000).
   DOI: [10.1006/jcph.2000.6575](https://doi.org/10.1006/jcph.2000.6575)

**流体/MHD AMR**
7. A. L. Zachary, A. Malagoli, P. Colella, *A higher-order Godunov method for multidimensional
   ideal magnetohydrodynamics*, SIAM J. Sci. Comput. 15(2), 263–284 (1994).
   DOI: [10.1137/0915019](https://doi.org/10.1137/0915019)
8. D. S. Balsara, *Divergence-free adaptive mesh refinement for magnetohydrodynamics*,
   J. Comput. Phys. 174(2), 614–648 (2001).
   DOI: [10.1006/jcph.2001.6917](https://doi.org/10.1006/jcph.2001.6917)
9. R. Teyssier, *Cosmological hydrodynamics with adaptive mesh refinement*,
   Astron. Astrophys. 385, 337–364 (2002).
   DOI: [10.1051/0004-6361:20011817](https://doi.org/10.1051/0004-6361:20011817)
10. J. M. Stone *et al.*, *The Athena++ adaptive mesh refinement framework: design and
    magnetohydrodynamic solvers*, Astrophys. J. Suppl. Ser. 249, 4 (2020).
    DOI: [10.3847/1538-4365/ab929b](https://doi.org/10.3847/1538-4365/ab929b)

**EB / cut-cell**
11. R. B. Pember, J. B. Bell, P. Colella, W. Y. Crutchfield, M. L. Welcome,
    *An adaptive Cartesian grid method for unsteady compressible flow in irregular regions*,
    J. Comput. Phys. 120(2), 278–304 (1995).
    DOI: [10.1006/jcph.1995.1165](https://doi.org/10.1006/jcph.1995.1165)
12. S. Osher, J. A. Sethian, *Fronts propagating with curvature-dependent speed: algorithms
    based on Hamilton–Jacobi formulations*, J. Comput. Phys. 79(1), 12–49 (1988).
    DOI: [10.1016/0021-9991(88)90002-2](https://doi.org/10.1016/0021-9991(88)90002-2)
13. S. Dey, R. Mittra, *A locally conformal finite-difference time-domain (FDTD) algorithm for
    modeling three-dimensional perfectly conducting objects*, IEEE Microw. Guided Wave Lett.
    7(9), 273–275 (1997). DOI: [10.1109/75.622536](https://doi.org/10.1109/75.622536)

**材料亚网格**
14. A. Farjadpour *et al.*, *Improving accuracy by subpixel smoothing in the finite-difference
    time domain*, Opt. Lett. 31(20), 2972–2974 (2006).
    DOI: [10.1364/OL.31.002972](https://doi.org/10.1364/OL.31.002972)
15. A. Deinega, I. Valuev, *Subpixel smoothing for conductive and dispersive media in the
    finite-difference time-domain method*, Opt. Lett. 32(23), 3429–3431 (2007).
    DOI: [10.1364/OL.32.003429](https://doi.org/10.1364/OL.32.003429)
16. T. Xiao, Q. H. Liu, *Enlarged cells for the conformal FDTD method to avoid the time step
    reduction*, IEEE Microw. Wirel. Components Lett. 14(12), 551–553 (2004).
    DOI: [10.1109/LMWC.2004.837384](https://doi.org/10.1109/LMWC.2004.837384)

**PIC / 电荷守恒**
17. J.-L. Vay *et al.*, *Mesh refinement for particle-in-cell plasma simulations: applications
    to and benefits for heavy ion fusion*, Laser Part. Beams 20(4), 569–575 (2002).
    DOI: [10.1017/S0263034602204139](https://doi.org/10.1017/S0263034602204139)
18. J.-L. Vay *et al.*, *Modeling of a chain of three plasma accelerator stages with the WarpX
    electromagnetic PIC code on GPUs*, Phys. Plasmas 28, 023105 (2021).
    DOI: [10.1063/5.0028512](https://doi.org/10.1063/5.0028512)
19. T. Zh. Esirkepov, *Exact charge conservation scheme for particle-in-cell simulation with an
    arbitrary form-factor*, Comput. Phys. Commun. 135(2), 144–153 (2001).
    DOI: [10.1016/S0010-4655(00)00228-9](https://doi.org/10.1016/S0010-4655(00)00228-9)

**矩/动理学**
20. J. Juno, A. Hakim, J. TenBarge, E. Shi, W. Dorland, *Discontinuous Galerkin algorithms for
    fully kinetic plasmas*, J. Comput. Phys. 353, 110–147 (2018).
    DOI: [10.1016/j.jcp.2017.10.009](https://doi.org/10.1016/j.jcp.2017.10.009)
21. L. Wang, A. H. Hakim, A. Bhattacharjee, K. Germaschewski, *Comparison of multi-fluid moment
    models with particle-in-cell simulations of collisionless magnetic reconnection*,
    Phys. Plasmas 22, 012108 (2015).
    DOI: [10.1063/1.4906063](https://doi.org/10.1063/1.4906063)

**教科书**
22. T. Plewa, T. Linde, V. G. Weirs (eds.), *Adaptive Mesh Refinement — Theory and Applications*,
    Springer (2005). DOI: [10.1007/3-540-27039-6](https://doi.org/10.1007/3-540-27039-6)
23. R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge University Press
    (2002)（重构、Riemann 求解、守恒的教科书标准参考）。
24. C. K. Birdsall, A. B. Langdon, *Plasma Physics via Computer Simulation*（PIC 噪声、色散与
    守恒；AMR 下的粒子噪声需按此书估计）。

---

## 7. 与配套验证文档的接口

- 本文的 R1–R8 对应配套文档 `benchmarks/EM_EB_MHD_VALIDATION_10CASES_v2.md` 的 §8.3：
  - P1（M2，介质 EB）→ 本文 R2/R5；
  - P2（M4，AMR 材料）→ 本文 R1/R4/R7；
- C10（核壳球 + AMR）是 M4 的验收算例；C2（斜平面介质 cut-cell）是 M2 的验收算例。
- 实施时**先 M2 后 M4**：M2 解决"介质 + EB 单层"，M4 才谈"每层"。

---

## 附录　文献检索与核对记录

- 所有 DOI 通过 `https://api.crossref.org/works?query.bibliographic=...` 检索并核对标题、期刊、年份；
- 书籍（22–24）核对版次与出版社；
- 代码事实（`FillPatcher` 存在但未用、`lev==maxLevel()` 限制、`build_eb_data_per_level`）
  均在 `<repo>` 与 `~/amrex-development` 中实测，
  行号见正文；
- 本文不含仿真数值。
