# Cut-cell 实现调研报告（项目 M2a 前期调研 · 不实现）

日期：2026-10-31
目标：为「完整 Maxwell + cut cell + 电介质 ε/μ/σ + 多流体」求解器厘清 cut-cell 的实现方法，确定 WarpX 已有/缺失，并给出选型问题。

---

## 1. 结论摘要：WarpX 已有 vs 缺失

| 模块 | WarpX 现状 | 评价 |
|---|---|---|
| B 场 cut-cell | **已有**：ECT（Enlarged Cell Technique，面扩展/面积借用）`EvolveBCartesianECT` | 真共形 cut-cell |
| E 场 cut-cell | **缺失**：仅 `eb_update_E` 掩膜做阶梯近似（跳过被覆盖边），无任何共形更新 | **核心缺口①** |
| 宏观 ε/μ/σ | **与 EB 正交**：普通 cell-centered `MultiFab`，不感知 cut-cell；无 ECT 分支 | **核心缺口②** |
| 流体 FV cut-cell | **缺失**：F1–F3 的流体求解器在规则网格上，未接 EB | **核心缺口③** |
| 小胞元稳定化 | B 场用 ECT 面积扩展 + BCK 兜底；**无 dt 缩减** | 几何侧稳定 |
| 几何原语 | 全部来自 AMReX `EBFArrayBoxFactory`（edge/face/volume fraction、质心、法向、边界面积） | 可用 |

---

## 2. Maxwell cut-cell（场求解器）

### 2.1 方法谱系
1. **Dey–Mittra 1997**（IEEE MTT-L 7:273–275）：局部共形 FDTD 起点——cut 边上用实际边长/面积做 Faraday/Ampere 的积分离散，避免阶梯误差。
2. **ECT**（Xiao & Liu 2008，IEEE TAP 56:765–773；Benkler–Chavannes–Kuster 2005 APS）：**Enlarged Cell Technique**——小 cut 面从稳定邻居「借面积」，把有效面积抬到稳定阈值 `S_stab`，消除 `1/S` 奇点导致的小胞元时间步缩减。WarpX 已实现（B 场）。
3. **电介质共形 FDTD**：Dey–Mittra 介质推广；Meep 亚像素有效 ε 张量（法向调和/切向算术）；Wang–Yin–Liu 2009（高阶 (2,4) 介质共形）。

### 2.2 WarpX 现有 ECT 细节（B 场，可复用）
- 几何量：`edge_lengths`（切向边长）、`face_areas`（cut 面面积）、`area_mod`（扩展后面积）、`borrowing`（借用关系）、`flag_info_cell`（extended/available/intruded/bck_stabilized）。
- `ComputeSStab`：稳定面积阈值 `0.5·max{邻接边长×垂直间距}`。
- `ComputeFaceExtensions`：单面借用 → 八邻域按面积分摊 → BCK 兜底（直接扩面积）。
- `EvolveECTRho`：Stokes 面环量 `∮E·dl/S`（注意：不是电荷重分配，名字误导）。
- `EvolveBCartesianECT`：`B -= dt·Rho_enl`，扩展胞元总环量在扩展面积上归一化。

### 2.3 E 场共形更新（缺口，需新写）
E 场目前是阶梯掩膜，需补 **Ampere 定律的 cut-cell 对偶离散**：用 `face_areas`/`edge_lengths`/`face_cent` 对 `∂E/∂t = c²∇×B − J/ε0` 做共形积分更新，并套用与 B 场对称的 ECT 面积扩展。这是 M2a 最大的新代码。

---

## 3. 流体 FV cut-cell（小胞元稳定化）

### 3.1 问题
显式 FV 的 CFL：`dt ≤ κ·dx/λ`，κ=体积分数。cut cell 的 κ→0 使 dt 无界缩小。

### 3.2 方法（三类，均有成熟实现）
1. **h-box**：Helzel–Berger–LeVeque，SIAM J. Sci. Comput. 26 (2005) 785–809；Berger–Helzel–LeVeque "h-box methods"（把 cut cell 与邻居合并成「h-box」做更新）。
2. **通量重分配（flux redistribution）**：小胞元超限通量按权重分给邻居（Colella/Pember 系）。
3. **状态重分配（state redistribution, SRD）**：Berger & Giuliani，JCP 2020；weighted SRD，SISC 2024。守恒、保线性、不降阶、显式、按全胞元体积取 dt。

### 3.3 AMReX 已自带（重大复用资产）
- `AMReX_EB_Redistribution.{H,cpp}`：`single_level_redistribute` / `single_level_weighted_redistribute` / `apply_flux_redistribution`（flux redistribution，含粗-细）。
- `AMReX_EB_StateRedistribute.cpp` + `AMReX_EB_StateRedistUtils.cpp` + `...SlopeLimiter_K.H`：state redistribution。
- `AMReX_EBFluxRegister.{H,cpp}`：EB 感知的 AMR reflux（M4 用）。

→ 流体小胞元稳定化**几乎不用自写**，可直接复用 AMReX 的 SRD 或 flux redistribution。

---

## 4. 电介质 ε/μ/σ 在 cut-cell（缺口，需新写）

现状：`MacroscopicProperties` 是普通 cell-centered MultiFab，不感知 EB。
可选离散：
1. **共形介质 FDTD（Dey–Mittra 介质版）**：cut 面用实际面积/边长加权 ε，切向 E、法向 D 连续。
2. **有效 ε 张量（Meep 亚像素）**：法向调和、切向算术，界面体元构造各向异性 ε⁻¹。
3. **体积分数加权 ε + 界面嵌入 Riemann**（与之前定的「真 cut-cell 界面条件」一致）。
σ 用 CN 隐式 `1/(1+σdt/2)`；μ 同理。

---

## 5. 文献清单

| 主题 | 文献 |
|---|---|
| 共形 FDTD（PEC） | Dey & Mittra, IEEE Microwave Guided Wave Lett. 7:273–275 (1997) |
| ECT | Xiao & Liu, IEEE Trans. Antennas Propag. 56:765–773 (2008)；Benkler–Chavannes–Kuster (2005) |
| 介质共形 | Wang, Yin, Liu, Electronics Letters 45(3) (2009)；Meep subpixel smoothing（Farjadpour 2006/2009, Kottke） |
| 小胞元 FV | h-box：Helzel–Berger–LeVeque, SIAM JSC 26:785–809 (2005) |
| 状态重分配 | Berger & Giuliani, JCP (2020)；weighted SRD, SISC 46:A2848 (2024) |
| 通量重分配 | Colella/Pember 系（embedded boundary Cartesian grid） |
| AMReX 实现 | `AMReX_EB_Redistribution`、`AMReX_EB_StateRedistribute`、`AMReX_EBFluxRegister` |

---

## 6. 关键实现路径（初步）

- **场**：复用 WarpX 的 B 场 ECT；新写 E 场共形 Ampere 更新（`face_areas`/`edge_lengths`/`face_cent` + ECT 面积扩展）；ε/μ/σ 改为 EB-aware（面/体积加权）。
- **流体**：F1–F3 的 FV 驱动接 EB 体积分数，小胞元用 AMReX SRD（或 flux redistribution）稳定。
- **AMR（M4）**：复用 `AMReX_EBFluxRegister` + WarpX 的 EB 感知 AMR。

---

## 7. 最终决策（已确认 2026-10-31）

| 决策点 | 定稿 |
|---|---|
| E 场 cut-cell | **Dey–Mittra 共形 Ampere 更新 + ECT 面积扩展稳定**（与 B 场对称，完整共形） |
| B 场 cut-cell | 复用 WarpX 已有 ECT（`EvolveBCartesianECT`） |
| 流体小胞元 | **复用 AMReX `AMReX_EB_StateRedistribute`（state redistribution, SRD）** |
| 电介质 ε/μ/σ | **M2a 先做 PEC，介质留到 M2b** |
| M2a 交付 | **PEC 共形闭环**：E+B cut-cell + 流体 SRD，验证绕射/反射 |

---

## 8. M2a 设计规范（B+C 合并版 · 2026-10-31 定稿）

> 依据：12 问答复（2D/3D 生产级 TE+TM、合并设计、可改 fork、扩展/借用、无 RZ、覆盖 AMR、
> 2D TE 与 3D 逐点一致）+ `benchmarks/CUTCell_LITERATURE_RESEARCH.md`（Xiao–Liu 2008 成对构造、
> Dey–Mittra 1997、Berger–Giuliani 2021 SRD、issue #7206 核对）。

### 8.0 合并的核心认识

B（2D TE）与 C（E 场共形 Ampere）**是同一个"对偶几何"问题的两个侧面**：
- 3D ECT 的 B 侧用 B 面面积 + E 边长；E 侧（C）用 **E 边的对偶面面积 + B 边的对偶边长**。
- 2D XZ 下：TM 的 B 侧 = 现有 ECT（已验证）；**TE 的 B 侧 = 3D x/z 面的 per-unit-y 约化**；
  E 侧（C）在 2D 同时覆盖 TM（Ex/Ez 对偶面）与 **TE（Ey 节点对偶面积）**。
- 因此 E/B 两侧、TM/TE、2D/3D 共用一套"对偶几何"基建，一次做完。

### 8.1 目标
在 fork（`~/WarpX_fluid`，**允许修改**）实现 **PEC 上 E+B 全分量共形 cut-cell Maxwell（2D XZ + 3D，TM/TE 常开、不切换）+ 流体 FV + SRD**。介质归 M2b。

### 8.2 场侧（B+C 合并）

**(a) 2D per-unit-y 几何（B 的前置）**
- `face_areas[0]` = x-法向线段开放长度（`getAreaFrac()[0]×dz`，即 Bx 位置）；`face_areas[2]` = `getAreaFrac()[1]×dx`（Bz 位置）
- `edge_lengths[1]` = 节点开放度（1/0）
- `EvolveB.cpp:261` ECT 循环 `idim<AMREX_SPACEDIM` 扩为 `idim<3`（否则 Bz 永不更新）
- `EvolveRhoCartesianECT` 补 XZ 的 `Rhox`/`Rhoz`（Ey 项，Ez 项 y-不变相消）
- **⚠ 范围修订（源码核查 §4）**：当前代码把上述量硬置零（`EmbeddedBoundaryInit.cpp:465/:535`），且 2D TE 实际需触碰 **12+2 处**（含 `ComputeSStab`/`MarkExtensionCells` dim=0 的 `k+1` 越界、`GetNeigh/SetNeigh` dim=0 的 3D-y 索引错位、BCK 的 `#if` 禁用、`shrink_borrowing` 与 `CountExtFaces` 归约只到 SPACEDIM）——完整清单见 `benchmarks/M2A_SOURCE_REVIEW_FINDINGS.md` §4

**(b) E 侧对偶几何共形 Ampere（C）**
- `EvolveECartesianECT`：`E += dt·(c²∮B·dl_dual − J·A_dual/ε0)/A_dual`
- 对偶面面积 ≠ 现有 `face_areas`（那是 B 面）；对偶量在**半格偏移网格**上（格心/边心），须新建数组（如 `dual_edge_length[3]`、`dual_face_area`）；对偶 `S_stab`/borrowing/`Venl` 式记账/BCK 兜底，按 Xiao–Liu 2008 成对构造
- **2D 分量语义（源码核查 §5.2）**：Ex/Ez 的切割体现在**对偶边长**（对偶面 y 向无限、A_dual=1/每单位 y）；Ey 的切割体现在**对偶面积**——节点对偶 cell 为 4 邻胞各贡献 **¼** 开放面积，即 `A_dual = Σ(¼)·A_open`（**非** 4 个面积直接相加，否则大 4 倍）
- J 项（源码核查 §5.3）：现 `EvolveE` 直接减 `μ0·j`；一期先做 **J=0**（真空+PEC 闭环），对偶面 J 积分与 ∇·E−ρ 一致性随 M2a-2 流体接入再验
- 掩膜语义保留（covered 区 E_t=0），须与对偶面积为零的处理一致

**(c) 稳定化**：两侧统一用 **扩展/借用（ECT 式）**（答复 8）；BCK 兜底前先核对上游 issue #7206

### 8.3 流体侧
- F1–F3 的 FV 驱动接 EB（`vol_frac` 加权通量）+ AMReX SRD；源项（Lorentz）**算子分裂**单独验收守恒
- PEC 壁面反射：ghost-cell + 镜像 Riemann（Colella 2006 做法）

### 8.4 验证（吸收文献报告 §5）
1. 2D TM 基线 r=−0.9993（**改动前基线**：换全共形 E 后该值会变、预期更接近 −1；新旧分开记录，不作回归目标）
2. **2D TE 与 3D 逐点一致**（答复 11，硬判据）：3D 对照设置 = **1 层 y + y 向周期 + 同 dx/dz/dt**（dy 与 per-unit-y 约定对齐）；容差建议逐点相对 ≤1%（待确认）
3. 3D 斜 PEC s/p 偏振：|r+1|≤1%、相位≤2°、伪散射≤1e-3
4. PEC 球 RCS vs Mie ≤2%（ka≥1）；2D 柱（TM）≤3%；腔模式 ≤0.5%（复用 rotated_cube）
5. 收敛阶 ≥1.5（3 分辨率）；能量 100 周期漂移 ≤1e-3
6. 流体：SRD 无源守恒 ≤1e-10，带源 ≤1e-6；流体冲击 PEC 壁
7. 单测：规则网格上"共形格式退化回 Yee"逐位一致；"借出通量=借入通量"

> **L2 判据订正（2026-10-31，见 `ECT_2D_TE_REPORT.md` §14）**：
> "规则网格"必须理解为**无任何 cut cell**（例如把 EB 面移出域外），已实测 201 步逐位一致。
> 与网格**对齐**的 PEC 面**不算**规则网格：AMReX 会把贴面单元标成 cut（`volfrac=1` 但 `isRegular()=false`），
> WarpX 的 staircase 掩膜据此冻结一整层，而共形几何按开放处理 ⇒ 两者等效壁面相差一层 `dz`。
> 因此"对齐面逐位等于 Yee"不可作为判据；替代判据为 L2a（规则域逐位）+ L2b（等效壁面位置 vs 解析）。

### 8.5 实现子任务
- **M2a-1a 对偶几何与判据**（2D per-unit-y + E 侧对偶面/边 + 对偶 `S_stab`；先做"退化回 Yee"单测与斜平面解析对偶面积对比）
- **M2a-1b 更新内核**（Rhox/Rhoz XZ 分支 + 循环扩 idim<3 + `EvolveECartesianECT` + 记账 + BCK）
- **M2a-2 流体**（EB 体积分数 + SRD + 源项分裂）
- **M2a-3 验证**（§8.4 全套）

### 8.6 约束落实
可改 fork（答复 4）；TM/TE 常开不切换（5）；场位置不动、保持标准 Yee（6/7）；扩展/借用（8）；不做 RZ（9）；**设计覆盖 AMR、先单层实现**（10，多层归 M4 用 `AMReX_EBFluxRegister`）。

> 下一里程碑 M2b = 电介质（subpixel 先行 → 共形介质对照，见文献报告 D4），M4 = AMR。
