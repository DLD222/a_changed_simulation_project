# C 侧调研：E 场对偶共形 Ampere 的文献依据、可修复性与实施框架

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 日期：2026-10-31 · 状态：**调研完成，待实现**
> 目的：回答三个问题——(1) 文献上 C 侧该怎么做；(2) C 能否修复已检出的问题；
> (3) 实施框架与步骤流程是什么。
> 前置：`C_SIDE_DESIGN.md`（公式细节）、`M2A_CODE_AUDIT.md`（问题清单）、`HANDOFF.md`（当前判据）。

---

## 0. 摘要（先给结论）

**C 侧 = 把 Ampere 定律也写成"对偶面环量"**，与已完成的 B 侧（Faraday，`EvolveBCartesianECT`）对称。
文献上这不是新方法：它就是 **Dey–Mittra 局部共形 FDTD 中"E 更新用开放面积"的那一半**，WarpX 的
ECT（Enlarged Cell Technique）只实现了 Faraday 那一半。

**能不能修复？** 分问题看，不能一概而论：

| 已检出的问题 | C 能否修复 | 判据/依据 |
|---|---|---|
| 米氏收敛阶测不出（`dx≤0.1` 数十步爆到 1e306） | **很可能能**（必要条件，非充分条件） | Dey–Mittra 双侧共形 + BCK 型稳定化正是为该问题设计；Nieter 2009 在 VORPAL 上验证过 Dey–Mittra 的 3D 可用性 |
| 45° 斜面 TM 长跑 1.9% 能量增长 | **预期能** | 单侧共形破坏离散能量的伴随性（§3.1），C 侧恢复对称 |
| L3 30° 掠入射 FAIL（斜率 −10.5%，TM 横向性 0.446） | **部分改善** | E 侧阶梯的 O(dx) 误差被去掉；但"反射/入射脉冲在探针处重叠"是判据缺陷，需同时改判据 |
| L2b 等效壁面偏 +0.37 ns（约差一层 `dz`） | **不能** | 根因是 AMReX 把"贴面单元"标为 cut（`volfrac=1` 但 `isRegular()=false`），属几何标记语义，与 E 侧精度无关；且 C 侧上线后该基线**会变**，必须重测 |
| AMR 多层几何未初始化 | **不能** | `InitializeEBGridData` 的 `lev == maxLevel()` 限制，是独立的初始化流程问题（§6 of `C_SIDE_DESIGN.md`） |
| 2D/3D 逐点一致性未知 | **不能** | per-unit-y 约化与 3D 面环量在 cut cell 上无先验等价性 |
| 将来加粒子后的电荷守恒 | **不能自动解决** | cut cell 上的电流沉积需单独的守恒格式（§2.5） |

**一句话**：C 是"把 1 阶阶梯变成可收敛格式"的**必要条件**，也是米氏爆炸与 TM 能量增长的**唯一根因级修复**；
但它不解决几何标记、AMR、PIC 电流沉积这三类问题，不应期待"做完 C 一切变好"。

---

## 1. 要修复的问题（来自既有审计，带原始数字）

| 编号 | 问题 | 证据 | 归属 |
|---|---|---|---|
| P-C1 | E 更新仍是掩膜 Yee：`m_eb_update_E` 只决定"是否更新"，不改变差分权重 | `EvolveE.cpp:90-98` 中 `ECT` 与 `Yee` 走同一个 `EvolveECartesian<CartesianYeeAlgorithm>`；`grep EvolveECartesianECT` 无命中 | 本文件主题 |
| P-C2 | 单侧共形导致小胞元局部增益 | 45° 斜面 TM 10000 步 `Emax/E0 = 1.019`，TE 仅 `1.0002` | §3 |
| P-C3 | 曲边界极小胞元导致数值爆炸 | 圆柱 `dx ≤ 0.1` 时 `volfrac ~ 1e-3`，80²/160² 在数十步内到 `1e306–1e308`；CFL 0.4/0.2/0.1 均爆 | §4.2 |
| P-C4 | 斜入射精度受限于 E 侧阶梯 | L3：45°/60° PASS（平均 `|r∓1| = 0.026–0.051`），30° FAIL（斜率 10.5%，TM 横向性 0.446） | §4.3 |
| P-C5 | 能量守恒无判据 | 尚无"真空腔 100 周期漂移 ≤ 1e-10"类测试（`C_SIDE_DESIGN.md` 的 C-5 未执行） | §5.6 |

---

## 2. 文献调研

### 2.1 理论骨架：Yee 格式是离散 de Rham 复形上的 `d` 与 Hodge star

把电场看作**边 1-形式**、磁场看作**面 2-形式**，则半离散 Maxwell 为

```
B^{n+1/2} = B^{n-1/2} − dt · (⋆₂)^{-1} C E^n          (Faraday，C 是离散 curl)
E^{n+1}   = E^n     + dt · (⋆₁)^{-1} C^T B^{n+1/2}     (Ampere，C^T 是对偶 curl)
```

`C` 只含 ±1（拓扑量，与几何无关），所有几何信息都在两个 Hodge star
`⋆₁`（边→对偶面）与 `⋆₂`（面→对偶边）里。**规则网格上两者互为转置的逆**：
`⋆₂ = diag(A_face)`、`⋆₁ = diag(A_dual)`，且 `A_dual · A_face` 在每个 cell 上是常数 `V²`，
因此 `C^T = −C`（在适当内积下）给出严格的离散能量守恒。

cut cell 上，一旦你**只改 `⋆₂`（面面积）而不改 `⋆₁`（对偶面积）**，`C` 与 `C^T` 的伴随性被破坏，
离散 Poynting 定理多出一项 `E·(⋆₁^{-1} − ⋆₂) C E`，其符号不固定 → 要么增长要么耗散。
这就是"单侧共形"在数学上的全部问题，与实现细节无关。

**教科书依据**：
- Taflove & Hagness, *Computational Electrodynamics*, 3rd ed., 2005 —— 第 3 章 Yee 格式、
  第 4 章 Von Neumann 稳定性分析（本项目的 `S_stab`/`A_stab` 判据都源自这一框架）；
- Bossavit, *Computational Electromagnetism*, 1998 —— 离散 Whitney 形式、对偶复形与 Hodge 算子，
  是"为什么 `C` 与 `⋆` 必须成对"的严格出处；
- LeVeque, *Finite Volume Methods for Hyperbolic Problems*, 2002 —— 守恒型离散与
  cut-cell 的通量形式（把 `A_dual` 理解为"守恒量的控制体"）。

### 2.2 Dey–Mittra（1997）：双侧共形是原始设计，不是可选优化

Dey & Mittra 的局部共形算法对**两侧同时**修正：

- **Faraday（H 更新）**：用被 PEC 截断后的**开放面积** `A_f` 作分母、**开放边长** `l_e` 作权重
  （PEC 上 `E_t = 0`，落在导体上的那部分环路贡献为零）；
- **Ampere（E 更新）**：用 E 边的**对偶开放面积** `A_dual`（围绕该边的四个 cell 的开放面积各取 1/4）。

即：`A_dual` 与 `A_f` 是同一次几何剖分的两侧，必须同时替换。WarpX 当前只做了前者。

**小胞元问题**：Dey–Mittra 原文即指出，当 `A_f/A` 极小时，显式格式的稳定步长 `dt ∝ √(A_f)` 会
任意缩小；工程做法是引入阈值 `α`（常取 1e-2 ~ 5e-2），把 `A_f/A < α` 的 cell 直接"丢弃"
（其 E 点强制为零或从更新中移除），以牺牲局部精度换取全局稳定。

**验证性实现**：Nieter, Cary, Werner, Smithe, Stoltz, *J. Comput. Phys.* **228** (2009) 7902–7916
把 Dey–Mittra 移植到 VORPAL，明确指出"最大稳定时间步可以任意小"是该方法的固有代价，
必须配合阈值/稳定化。这是本项目 P-C3 的直接文献对应。

### 2.3 稳定化的三条路线（WarpX 已在用第 2 条的一半）

| 路线 | 做法 | 代价 | 代表 |
|---|---|---|---|
| 丢弃小胞元（Dey–Mittra 阈值） | `A_f/A < α` 的 cell 的 E 置零 | 局部精度降到 0 阶；α 敏感 | Dey–Mittra 1997；Nieter 2009 |
| **BCK** | 把 `A_f` 抬到 `A_stab`（`S_stab`），使之不缩小稳定步长 | 精度略降；需给出 `A_stab` 的闭式稳定性判据 | Benkler, Chavannes, Kuster, *IEEE TAP* **54**(6) 1843–1849, 2006, DOI 10.1109/TAP.2006.875909 |
| 扩胞 / USC | 把不规则 cell 扩大到邻居，或构造一致稳定的共形格式（USC） | 实现复杂；需邻接记账 | Zagorodnov, Schuhmann, Weiland, *J. Comput. Phys.* **225** (2007) 1493–1507；ECT = Enlarged Cell Technique（WarpX glossary, IEEE 4463918） |

**关键推论**：WarpX 的 ECT 名字就来自"扩胞"。扩胞的本质是**人为改变几何量以换取稳定性**。
如果只在 Faraday 侧扩（`S → S_stab`）、Ampere 侧不扩（`A_dual` 还是"真实"值），
那么两侧描述的**就不再是同一块几何** —— 这正是 P-C2 的机理。C 侧必须成对扩。

### 2.4 开源项目对照

| 项目 | 共形方式 | 对 PEC cut cell 的处理 | 对本项目的意义 |
|---|---|---|---|
| **WarpX / AMReX EB** | `face_areas` + `area_mod` + `ECTRhofield` | 只有 Faraday 侧共形，Ampere 侧是掩膜 Yee | 本项目要补的就是这一半；上游无可参照实现 |
| **Meep** | subpixel smoothing（材料界面按 `ε`/`1/ε` 平均、各向异性张量） | PEC 边界不是 subpixel 的目标场景；文档明确"直接离散不连续界面会把 2 阶降为 1 阶" | 提供"阶梯化 ⇒ 1 阶"的权威陈述与开源证据；但其平均是对**材料系数**，不能直接搬给 PEC 的**几何面积** |
| **VORPAL / VSIM** | 完整 Dey–Mittra（双侧 + 阈值） | 直接实现双侧共形 | 唯一有公开论文的"双侧共形 PIC"参考实现（Nieter 2009） |

Meep 的 subpixel 平滑（Farjadpour et al., *Opt. Lett.* **31** (2006) 2972）与本文问题同源
（都是"界面与网格错位导致 O(dx) 误差"），但它改的是**材料张量**；PEC 的共形必须改**几何面积/边长**，
因为 PEC 的物理约束是"切向 E = 0"这一**边界条件**，而非系数跳变。这是两条不能互换的技术路线。

### 2.5 PIC 特有：电流沉积（本阶段不涉及，但要预先知道）

本次米氏测试无粒子（`J = 0`），C 侧只需处理 `∂E/∂t = c²∇×B`。
一旦加入粒子，cut cell 上的电荷守恒要求电流沉积沿**同一条开放边**记账：
`Δt` 内穿过开放面的通量必须与 `ρ` 的开放体积变化严格抵消。
Dey–Mittra/VORPAL 与后续的共形电荷守恒沉积（见 IEEE 的 Conformal EM-PIC 综述，Xplore doc 7302583）
都需为此额外改动 Esirkepov 型沉积。**这不在 M2a 范围内，但 C 侧的数据结构应预留"开放面积"给沉积使用**。

---

## 3. 理论分析：单侧共形为什么会失稳（可定量）

### 3.1 能量法（一步映射）

定义离散能量 `U = ½ Σ_c V_c |E|² + ½ Σ_f A_f |B|²`（`V_c` 为 cell 开放体积）。
规则网格上一步更新严格保持 `U`（除边界通量）。cut cell 上若 `⋆₁` 未共形：

```
ΔU/U ~ CFL · (A_dual^Yee / A_f − 1) · (局部 |E|² / 平均 |E|²)
```

- `A_dual^Yee / A_f > 1` ⇒ 每步注入能量（增长）；
- 比值 < 1 ⇒ 每步耗散（表现为"反射系数偏低"，即 L3 中 TE 64² 的 `r = −0.9188`）。

**这解释了一个此前未被解释的现象**：为什么 TE 在小胞元上"耗散"（`r` 偏向 0）而在极小球面胞元上"爆炸"。
两种表现的差别只在失配比的方向与幅度，机制相同。

### 3.2 toy 估计（量级对照，非精确预测）

取单个 cut cell，`f = A_f/A`。按 `⋆₁` 未共形的 Yee 更新，每步增长因子近似

```
g ≈ 1 + CFL · (1/f − 1)
```

| `f` | CFL | `g`（每步） | 从 `1e6` 到 `1e306` 所需步数 |
|---|---|---|---|
| 0.5 | 0.4 | 1.4 | ~2000 |
| 1e-1 | 0.4 | 4.6 | ~450 |
| 1e-2 | 0.4 | 40.6 | ~160 |
| 1e-3 | 0.4 | 400 | ~110 |

实测：80²（`dx = 0.1`，`volfrac` 最小约 `1e-3`）在**数十步**内爆炸，160² 更快。
toy 估计偏保守（实际还有 `ρ`/`A_stab` 的二次放大与多胞元叠加），但量级一致。
**结论：爆炸不是"CFL 不够小"，而是"两侧几何不一致"；降 CFL 无效（实测 CFL 0.1 仍爆）与推导一致。**

### 3.3 为什么 TE 只涨 1.9%、TM 涨得多

- TE：只有 `Ey` 一个分量参与，`A_dual = ¼ dx dz Σvolfrac`，失配是"面积比"；
- TM：`Ex`/`Ez` 用**对偶边长**（`L·dz`、`L·dx`）作分母，且 `Rhoy` 同时含 `Ex`、`Ez` 两项，
  两个分量的失配**相加**，故增长更快；
- 另外 TM 的 `Rhoy` 在 cut cell 上被 BCK 抬过（`Sy → S_stab`），而 `Ex`/`Ez` 的边长没有对应的抬升，
  这是"半侧稳定化"的又一处不一致。

**这给出 C 侧的一个额外判据**：C 侧实现后，TE 与 TM 的能量漂移应回到**同一量级**（差异 ≤ 一个数量级）。

---

## 4. C 侧的框架

### 4.1 算子定义（2D XZ，先做这一个）

沿用 `C_SIDE_DESIGN.md` §3.2/§4.1，此处只补"为什么是这四个量"：

| E 分量 | 对偶面（守恒控制体） | 对偶面积 `A_dual` | 对偶边长（乘在 B 上） |
|---|---|---|---|
| `Ex`（x 边） | 围绕该边、法向 `y` 的矩形（2D 下退化为沿 z 的线段） | `1`（每单位 y） | `Lz(i,j) = areafrac_z(i,j)·dx` |
| `Ez`（z 边） | 同上，法向 `y` | `1` | `Lx(i,j) = areafrac_x(i,j)·dz` |
| `Ey`（节点） | 围绕该节点的 4 个 cell 各 1/4 | `¼ dx dz Σ volfrac` | `Lx`、`Lz`（同 cell 面） |

`Ey` 的 `¼` 是最容易写错的一处：不是"4 个面积相加"，而是**每个 cell 只贡献 1/4**
（因为该 cell 只占节点对偶面的一个象限）。

### 4.2 数据结构与生命周期

```
新增（仅 ECT 且 dim==2 时分配，其余维度不分配，保证 C-0 逐位不变）：
  dual_face_area[lev]   : MultiFab(1 comp, IntVect(1,1))        节点居中
  dual_edge_length[lev] : 3 × MultiFab(1 comp, E_idim 交错)     边居中
生命周期与 face_areas 完全一致（同一处分配、同一处 FillBoundary、同一处释放）。
```

不在 `Fields.H` 里加新 `FieldType` 也可以（B 侧 `face_areas` 就在 register 里）；
但为了 AMR 的 `FillBoundary` 与 plotfile 调试，建议**加**，与 `face_areas` 并列。

### 4.3 开关与退化（最重要的一条工程纪律）

```
algo.ect_conformal_E = 0   (默认)
```

- 关闭时：代码路径与当前分支**逐位相同**（C-0）。
- 打开 + 规则域：必须与 Yee **逐位相同**（C-1，`max|d| = 0`）。
  实现顺序必须是：**先写死常数几何**（`A_dual ≡ dx·dz`、`L ≡ dz/dx`）跑通 C-1，
  再接真实几何量再跑一次 C-1。两次都要 0 才继续。

理由（来自 B 侧的经验）：索引半格、`area_frac` 的左右归属、2D 下 `k` 的哑元，
这三个坑任意一个都会让退化不成立，而"不成立"的表现是"看起来还行但差 1e-2"。

### 4.4 对偶稳定化

对称于 B 侧：`A_dual < A_stab` 时抬到 `A_stab`；`L < L_stab` 同理。

```
A_stab = 0.5 · max{ 相邻对偶面的 (L · Δ) }        (与 ComputeSStab 同构)
L_stab = 0.5 · min{ 邻接 cell 的完整边长 }
```

**不要**直接搬 B 侧的 `Venl` 借用记账：E 侧没有"环量累加器"的等价物，
先 BCK-only 做 `volfrac` 扫描（`1e-1 → 1e-6`），确认稳定性边界，再决定是否需要借用。

### 4.5 与现有代码的接缝（7 处，全部有行号）

| # | 文件 | 位置 | 动作 |
|---|---|---|---|
| 1 | `Source/Fields.H` | `FieldType` 列表（`:71` 附近） | 加 `dual_face_area`、`dual_edge_length` |
| 2 | `Source/WarpX.cpp` | ECT 分配段（与 `face_areas`/`area_mod` 同一处） | 按 §4.2 分配，仅 2D 且开关打开 |
| 3 | `Source/EmbeddedBoundary/EmbeddedBoundaryInit.{H,cpp}` | 新增 `ComputeDualGeometry()` | 由 `volfrac`/`area_frac` 填，规则域打印 min/max 自检 |
| 4 | `Source/Initialization/WarpXInitData.cpp` | `InitializeEBGridData`（ECT 分支） | 调用 `ComputeDualGeometry` |
| 5 | `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp` | `EvolveE`（`:55`）分派处 | `m_fdtd_algo==ECT && conformal_E` → 新内核 |
| 6 | 同上 | 新增 `EvolveECartesianECT` | 见 §4.6 |
| 7 | `Source/EmbeddedBoundary/WarpXFaceExtensions.cpp` | `ComputeSStab` 旁 | 新增对偶 `A_stab` 计算 |

**不改**：`EvolveB.cpp`（B 侧已正确）、`EvolveECTRho.cpp`、掩膜函数 `MarkUpdateECellsECT`
（共形后掩膜只应保留"covered 点强制 0"的语义，不再承担精度职责；这一语义调整放在 S6 之后）。

### 4.6 内核骨架（2D `Ey` 为例）

```text
A(i,j)·(Ey_new − Ey_old) =
    c² dt · [ (Bx(i,j)·Lz(i,j)   + Bx(i−1,j)·Lz(i−1,j))   / 2
            − (Bx(i,j−1)·Lz(i,j−1) + Bx(i−1,j−1)·Lz(i−1,j−1)) / 2
            − (Bz(i,j)·Lx(i,j)   + Bz(i,j−1)·Lx(i,j−1))   / 2
            + (Bz(i−1,j)·Lx(i−1,j) + Bz(i−1,j−1)·Lx(i−1,j−1)) / 2 ]
```

实现注意：
1. `A` 与 `L` 从 MultiFab 读，**不允许**在核心里现算 `dx*dz`（否则退化不可验证）；
2. `A` 有 1 层 ghost，`FillBoundary` 在 `ComputeDualGeometry` 后调用一次；
3. `A ≤ 0` 的节点跳过（covered），但**必须**在同一处把 `Ey = 0` 写回，否则会出现"零除保护留下脏值"；
4. `J` 项本阶段置 0（真空 + PEC），但保留接口形参，避免后续改签名。

---

## 5. 分步实施流程（S0–S10）

每步都有**通过判据**与**失败时怎么做**。没有通过判据不许进下一步（这条比任何公式都重要）。

| 步 | 内容 | 通过判据 | 失败时 |
|---|---|---|---|
| **S0** | 冻结基线：新分支 `feat/ect-conformal-E`；记录当前 8 组判据数字（L1–L4、长跑、回归） | `git status` 干净；基线表写入 `C_SIDE_RESEARCH.md` 附录 | — |
| **S1** | `ComputeDualGeometry`：只计算+打印，不接内核 | 规则域 `A_dual ≡ dx·dz`、`Lx ≡ dz`、`Lz ≡ dx`（`max|dev| < 1e-14`）；cut cell 上 `0 < A_dual ≤ dx·dz` | 查 `volfrac` 的 `IntVect` 是否节点居中、四象限索引是否 `(i−1|i, j−1|j)` |
| **S2** | 内核 + **常数几何**（`A≡dxdz`、`L≡d`），开关打开 | L2a：201 步 `max|d| = 0`（逐位） | 符号/左右归属错；用单步 `Ey/Bx` 手算比对 |
| **S3** | 接真实几何量 | L2a 仍 `max|d| = 0` | 索引半格错；打印 `A`/`L` 与解析值逐点比 |
| **S4** | L2b 重测（64/128/256） | 等效壁面误差 ≤ 0.5 ns 且不劣于当前（`+0.371/+0.017/+0.076`） | 允许基线变化，但必须**解释**变化方向 |
| **S5** | L3 45°/60° 复测；30° 复测（同时改判据：用窄脉冲或加大域） | 45°/60° 平均 `|r∓1| ≤ 0.03`（不劣于 0.026–0.051）；30° 斜率误差 ≤ 5% | 若变差，回到 S3 检查 `A` 的符号约定 |
| **S6** | **圆柱 C12 + 米氏 40/80/160** | 80²/160² 不再爆炸（`max|E|` 有限且 < 10× 入射峰值）；40² 结果与当前一致 | 若仍爆：查 `A_stab` 是否生效（打印最小 `A_dual`）；查掩膜是否仍冻结 cut cell 内的 E 点 |
| **S7** | 对偶 BCK（`A_stab`、`L_stab`）+ `volfrac` 扫描 `1e-1…1e-6` | 全扫描无爆炸；能量漂移 ≤ 1e-6/千步 | 扫描中定位临界 `f`，反推 `A_stab` 系数 |
| **S8** | TM 长跑 10000 步复测 | `|Eend−E0|/E0 ≤ 1e-4`（当前 9.5e-3）；TE 不劣于 3e-4 | 若 TM 仍涨，查 `Rhoy` 与对偶边长的交互（§3.3） |
| **S9** | 米氏 5 档收敛阶（TE/TM，L2 环 1.5R–2.5R） | 报告实测阶数（预期从 ~1 阶提升；不预设 2 阶） | 阶数 < 1：查 RCS 外推与近场采样是否引入额外误差 |
| **S10** | 3D 推广 + 上游回归 | `rotated_cube` 2D/3D 不变（`6.8e-10` / `≤1e-13`） | 3D 三个分量对称性；先只开 2D 开关，3D 走旧路径 |

**S0 之后可并行**的独立任务（不阻塞 C）：AMR per-level 初始化（`C_SIDE_DESIGN.md` §6）、
2D/3D 逐点对照脚本。

---

## 6. 风险与预案

| 风险 | 概率 | 影响 | 预案 |
|---|---|---|---|
| 索引/半格约定错，退化不成立 | **高** | S3 卡住 | 常数几何先行（S2）；打印 `A`/`L` 与解析值逐点比对；用单 cell 手算 |
| 符号/左右归属错 | 高 | L2a 差 O(1) | 用"上游/下、左/右各两段 /2"的对称写法，只可能整体符号错，一次修正即好 |
| `A_dual` 定义与掩膜语义冲突 | 中 | 小胞元行为反常 | 明确：共形管"权重"，掩膜只管"covered 强制 0"；两者同时生效时掩膜优先 |
| 稳定化不足，S6 仍爆 | 中 | 阻塞米氏 | S7 的 `volfrac` 扫描定位临界值；必要时退到 Dey–Mittra 阈值（丢弃 `f<1e-3`） |
| 性能下降（多 3 个 MultiFab + 每步读 6 个数组） | 中 | 米氏 640² 变慢 | 先测 160² 的每步耗时；必要时把 `A`/`L` 预乘进系数数组（`dt·c²/A`） |
| TM 基线（`r = −0.999361`）变化 | 高 | 需重记基线 | 预期变化；变化后重跑 L3/TM 全套并更新基线表 |
| AMR 下 `FillBoundary` 跨层不一致 | 中 | AMR 结果错 | 与 S10 解耦；AMR 单独一步，用"细层覆盖区 vs 单层细网格"比对（≤1e-8） |

---

## 7. 工作量

| 项 | 代码量 | 含调试 |
|---|---|---|
| S1 几何量 | ~150 行 | 0.5 天 |
| S2–S3 内核 + 退化 | ~130 行 | 1 天 |
| S4–S5 判据复测 | 0（复用脚本） | 0.5 天 |
| S6–S7 稳定化 + 扫描 | ~90 行 | 1–1.5 天 |
| S8–S9 长跑 + 收敛阶 | 0 | 0.5 天（640² × 5 档算力为主） |
| S10 3D + 回归 | ~250 行 | 1.5 天 |
| **合计** | **~620 行** | **5–6 天** |

风险预算：索引类问题占 50% 以上调试时间（B 侧的实测经验）。

---

## 8. 参考文献

**教科书**

1. A. Taflove, S. C. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain
   Method*, 3rd ed., Artech House, 2005.（Yee 格式、Von Neumann 稳定性）
2. A. Bossavit, *Computational Electromagnetism: Variational Formulations, Complementarity,
   Edge Elements*, Academic Press, 1998.（离散形式、对偶复形、Hodge 算子）
3. R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge Univ. Press, 2002.
   （守恒型离散、cut-cell 通量）
4. C. K. Birdsall, A. B. Langdon, *Plasma Physics via Computer Simulation*, Adam Hilger, 1991.
   （离散连续性方程、电荷守恒）

**论文**

5. S. Dey, R. Mittra, "A locally conformal finite-difference time-domain (FDTD) algorithm for
   modeling three-dimensional perfectly conducting objects," *IEEE Microwave and Guided Wave
   Letters*, vol. 7, no. 9, pp. 273–275, 1997.（双侧共形 + 小胞元阈值）
6. S. Benkler, N. Chavannes, N. Kuster, "A new 3-D conformal PEC FDTD scheme with user-defined
   geometric precision and derived stability criterion," *IEEE Trans. Antennas Propag.*,
   vol. 54, no. 6, pp. 1843–1849, 2006. DOI: 10.1109/TAP.2006.875909（BCK；WarpX 的 `S_stab` 来源）
7. I. Zagorodnov, R. Schuhmann, T. Weiland, "Conformal FDTD-methods to avoid time step reduction
   with and without cell enlargement," *J. Comput. Phys.*, vol. 225, no. 2, pp. 1493–1507, 2007.
   （USC / 扩胞；ECT 的理论出处）
8. C. Nieter, J. R. Cary, G. R. Werner, D. N. Smithe, P. H. Stoltz, "Application of Dey–Mittra
   conformal boundary algorithm to 3D electromagnetic modeling," *J. Comput. Phys.*,
   vol. 228, no. 21, pp. 7902–7916, 2009.（VORPAL 实现；小胞元导致 dt 任意缩小）
9. A. Farjadpour, D. Roundy, A. Rodriguez, M. Ibanescu, P. Bermel, J. D. Joannopoulos,
   S. G. Johnson, "Improving accuracy by subpixel smoothing in the finite-difference
   time-domain," *Optics Letters*, vol. 31, no. 20, pp. 2972–2974, 2006.（Meep 亚像素平均；
   材料界面 ≠ PEC 边界，但"阶梯 ⇒ 1 阶"的结论同源）

**开源项目 / 文档**

10. WarpX `Docs/source/glossary.rst`：`ECT = Enlarged Cell Technique`（IEEE 4463918）、
    `BCK = Benkler-Chavannes-Kuster`（IEEE 1638381）。
11. WarpX `Source/EmbeddedBoundary/WarpXFaceExtensions.cpp:173`、`EvolveB.cpp:376`：
    源码注释直接引用 BCK 作为小胞元稳定化。
12. Meep 文档 *Subpixel Smoothing*：明确"直接离散不连续界面会把二阶降为一阶"。
13. AMReX `EBLevelGrid` / `CutCell` 的 `volfrac`、`area_frac` 定义（本项目几何量来源）。

> 备注：检索中另见 IEEE 的 "Conformal Electromagnetic Particle in Cell: A Review"
> （Xplore 文档号 7302583），是共形 EM-PIC 的综述，但本次未核实其卷期，故不列入正式引用，
> 仅在 §2.5 作为线索。

---

## 附录 A：S0 基线表（实现前必须先填，实现后逐项对照）

| 判据 | 当前值 | C 侧后 | 变化说明 |
|---|---|---|---|
| L1 几何退化 | PASS | | |
| L2a 规则域 | `max|d| = 0`（201 步） | | 必须仍为 0 |
| L2b 往返误差 64/128/256 | `+0.371 / +0.017 / +0.076` ns | | 允许变，须解释 |
| L3 45° TE/TM 平均 `|r∓1|` | `0.0427 / 0.0263` | | 目标不劣化 |
| L3 60° TE/TM | `0.0434 / 0.0507` | | |
| L3 30° TE 斜率误差 | `−10.5%` | | 目标 ≤5% |
| L3 30° TM 横向性 | `0.446` | | 目标 ≤0.15 |
| L4 圆柱 镜面/阴影 TE | `5.35` | | |
| L4 圆柱 镜面/阴影 TM | `1.42` | | |
| 长跑 TE 漂移（1e4 步） | `3.0e-4` | | 不劣化 |
| 长跑 TM 漂移（1e4 步） | `9.5e-3` | | 目标 ≤1e-4 |
| 米氏 80² 稳定性 | 爆炸（`1e306`） | | 必须有限 |
| 米氏 160² 稳定性 | 爆炸 | | 必须有限 |
| 上游回归 2D/3D | `6.8e-10 / ≤1e-13` | | 不变 |
