# M2b 实施方案调研：电介质 σ/ε/μ（张量）+ 界面处理 + Fresnel R/T 验收

> 日期：2026-10-31 · 前置：`architecture-final.md`（决策 3/4、M2b 验收）、`C_SIDE_LS_ECT_PLAN.md`（场侧现状）
> 范围：M2b = 电介质 σ/ε/μ + 界面处理；本文给出**方法、复用件、验证阶梯、实施顺序**。
> 结论先行：
> 1. **场侧（FDTD）的 M2b = 共形等效 ε/μ 张量（subpixel 平均）+ σ 的 CN 隐式**；这正是 Yu–Mittra 共形介质、
>    Hwang–Cangellaris 有效介电常数、Meep subpixel（Kottke 微扰理论）这一条成熟路线。
> 2. **"界面嵌入 Riemann"不是 FDTD 的必需品**：切向 E/H、法向 D/B 连续在 Yee 网格上由"等效系数平均"
>    （切向算术 / 法向调和）**自动满足**；Riemann 通量是 FV/DG（流体、Gkeyll 风格 Maxwell）的做法，
>    应在流体侧（M2a-2 的 PEC 镜像 Riemann 同族）使用。本文 §4.7 给出这个定位及其依据。
> 3. WarpX 已有 `MacroscopicProperties`（σ/ε/μ 场 + 交错）与 `MacroscopicEvolveE.cpp`（`H=B/μ`、σ 的
>    alpha/beta 系数）可复用；新工作集中在"ε/μ 的**共形等效张量**与 ECT 对偶几何的接缝"。

---

## 1. 目标与验收（对齐 architecture-final.md 决策 3/4 与 M2b 验收）

| 项 | 内容 |
|---|---|
| 一期 | 常数 σ/ε/μ（分区常数介质），随 M2b-2 升为标量场、M2b-3 升为张量场、M2b-4 色散张量（ADE） |
| σ 处理 | CN 隐式 `1/(1+σdt/2)`（决策 3 已定；WarpX 已有对应实现） |
| 界面 | 真 cut-cell 上等效 ε/μ（切向算术、法向调和，各向异性张量）；D/B 连续性由平均自动满足 |
| 验收 | **EM 波介质平板 R/T vs Fresnel 解析**：垂直/斜入射、s/p、Brewster、临界角/全反射、σ 衰减；\|R\|/\|T\| ≤1%、相位 ≤2°、分辨率收敛；曲面介质与 Meep 对照 |

---

## 2. 文献调研

| 文献 | 要点 | 对本项目的意义 |
|---|---|---|
| Yu & Mittra 2001, "A conformal finite difference time domain technique for modeling curved dielectric surfaces"（[Google Scholar](https://scholar.google.ro/citations?view_op=view_citation&hl=ro&user=43RLFSYAAAAJ&citation_for_view=43RLFSYAAAAJ:u5HHmVD_uO8C)） | 共形 FDTD 处理**曲面介质**：把 cut cell 上的 ε 修正为共形面积/长度加权的等效 ε | M2b 的直接方法出处（与 Dey–Mittra PEC 对称的介质版） |
| Hwang & Cangellaris 2001, "Effective permittivities for second-order accurate FDTD equations at dielectric interfaces"（[Google Scholar](https://scholar.google.com.co/citations?view_op=view_citation&hl=en&user=j3SZmOcAAAAJ&citation_for_view=j3SZmOcAAAAJ:W7OEmFMy1HYC)） | 二阶精度的等效 ε 规则：界面**切向**用算术平均、**法向**用调和平均 | "平均规则"的严格出处 |
| Kottke, Farjadpour & Johnson 2008, "Perturbation theory for anisotropic dielectric interfaces…", PR E 77, 036611（[arXiv:0708.1031](https://ar5iv.labs.arxiv.org/html/0708.1031)） | 各向异性有效张量 ε̃⁻¹ = P⟨ε⁻¹⟩ + (I−P)⟨ε⟩⁻¹（P = n̂n̂ᵀ，n̂ 为界面法向）；并给出 **σ 的对应处理**与微扰误差分析 | Meep subpixel 的理论公式，M2b 的核心公式（含 σ 平均） |
| Farjadpour et al. 2006, Opt. Lett. 31:2972–2974（DOI 10.1364/OL.31.002972） | subpixel 平滑把"阶梯界面的一阶误差"恢复为二阶 | 动机与收敛性依据 |
| Deinega & Valuev 2007, Opt. Lett. 32:3429（DOI 10.1364/OL.32.003429） | 斜/对角界面的 subpixel 处理 | 斜界面（本项目的斜面算例）补充 |
| Taflove & Hagness 2005；TUM 博士论文（[mediatum PDF](https://mediatum.ub.tum.de/doc/1760138/1760138.pdf)） | 色散介质 ADE（Debye/Drude/Lorentz）与 CN-σ 更新模板 | M2b-4 的模板 |
| gprMax（本机 `<local-ref>/gprMax`） | ADE/递归卷积色散实现（Debye/Drude/Lorentz）、材料稳定性报告 | M2b-4 的开源参照 |

---

## 3. 开源对照与复用件

### 3.1 Meep（本机源码，`<local-ref>/meep/src/anisotropic_averaging.cpp`）

- 公式就在代码里（:201 注释）：`P*minveps + (I-P)*1/meps = P*(minveps−1/meps) + I*1/meps`；
- 法向 `n̂` 来自 **level set 梯度**（`normal_vector`，:100 `gradient(normal_vector(ft,v))`）；
- 平均是**体积积分**（`volume` 类），初始化一次性算好，不进每步内核；
- 局限（`Subpixel_Smoothing.md`）：只对**无损无频散介质**严格成立；金属/色散只取 ε∞；界面上场仍一阶。
- 对本项目：**直接照抄公式与"level set 梯度定 n̂ + 逐 cell 积分"的构造**，PEC（ε=∞ 侧）按 ε_max 近似或独立处理。

### 3.2 WarpX（fork 已存在，可直接复用）

| 件 | 位置 | 可复用内容 |
|---|---|---|
| `MacroscopicProperties` | `Source/FieldSolver/FiniteDifferenceSolver/MacroscopicProperties/` | σ/ε/μ 三个 MultiFab 的分配、初始化（解析/常数）、交错 IndexType（`getsigma_mf/getepsilon_mf/getmu_mf`、`*_IndexType`） |
| `MacroscopicEvolveE.cpp` | 同目录 | 宏观 Maxwell 的 E 更新模板：`H = B/μ`（`FieldAccessorMacroscopic` 插值）、σ/ε 插值到 E 位置、**alpha/beta 系数**（σ CN 隐式） |
| 开关 | `WarpX.cpp:452-454`（`m_em_solver_medium == MediumForEM::Macroscopic`） | 输入侧（`algo.em_solver_medium`）已接通 |

**注意**：上游 `MacroscopicEvolveE` 是**无 EB** 的宏观路径；M2b 的新工作是把它与 **ECT/对偶几何** 结合（见 §4.5），而不是重写。

### 3.3 gprMax

- 色散（Debye/Drude/Lorentz）用 ADE 场扩展 + 递归卷积；材料稳定性报告给出能量/稳定性分析框架。
- 用于 M2b-4 对照。

---

## 4. 方法设计

### 4.1 材料场（一期→二期）

```
一期（M2b-1）：分区常数 ε/μ/σ（复用 MacroscopicProperties 的常数初始化）
二期（M2b-2/3/4）：标量场 → 张量场（3×3 对称，存 cell 中心，按 Kottke 公式由 level set 构造）→ 色散张量（ADE 极点）
```

### 4.2 等效 ε/μ 张量（M2b 的核心公式）

对每个 cut cell（或对偶面/线）：

```
ε̃⁻¹ = P ⟨ε⁻¹⟩ + (I − P) ⟨ε⟩⁻¹       （Kottke 2008）
μ̃ 同构；P = n̂ n̂ᵀ；n̂ = ∇φ/|∇φ|（level set 梯度，规则区退化为体积分数平均）
⟨·⟩ = 沿对偶几何的开放部分积分平均（与 M2a 的 level-set 求交同一套工具）
```

2D 各分量落位（与 ECT 对偶几何对齐）：

| 量 | 界面关系 | 等效系数落位 |
|---|---|---|
| `Ex`（x 边，分母 Lz） | 法向连续 = D_x | ε 取**法向**（调和）于对偶段 |
| `Ez`（z 边，分母 Lx） | 法向连续 = D_z | 同上 |
| `Ey`（节点，分母 A_y） | 切向连续 = E_y | ε 取**切向**（算术）于对偶面 |
| `B` 侧 | 法向 B 连续 | μ 取法向（算术/调和对应侧） |

规则域退化为 Yee 的 `ε(i,j+½)`/`ε(i,j)` 平均；`L2a 逐位`判据沿用（μ=1 时与 PEC 版一致）。

### 4.3 σ 的 CN 隐式

- 更新形式：`E^{n+1} = α·E^n + β·(∇×H 项)`，`α = (1 − σdt/2ε)/(1 + σdt/2ε)`、`β = dt/(ε(1 + σdt/2ε))`；
  **`MacroscopicEvolveE.cpp:186-187` 的 `T_MacroAlgo::alpha/beta(sigma, eps, dt)` 已是这个实现**，直接搬到 ECT 内核。
- σ 的界面平均按 Kottke 2008 的 σ̃ 规则（与 ε̃ 同构：切向算术、法向调和），落位与 ε 相同。
- 稳定性：CN 对大 σ 无条件稳定，避免了"σ 显式要 dt < 2ε/σ"的限制。

### 4.4 μ 与 Faraday

- 保持 B 为基本量（ECT 的 `Rho` 不变）；Ampere 前算 `H = B/μ̃`（`FieldAccessorMacroscopic` 同款插值），
  E 更新用 H。Faraday 侧不变（μ 只进 H）。
- 界面：μ 的法向连续由平均自动满足（与 ε 对偶）。

### 4.5 与 ECT/对偶几何的接缝（真正的新代码）

| 处 | 改动 |
|---|---|
| `ComputeDualGeometry` 家族 | 新增 `ComputeEffectiveEpsilon/Mu`：沿对偶线/面对 level set 积分 ⟨ε⟩、⟨ε⁻¹⟩、n̂ → ε̃（张量或按分量预投影）；一次性初始化 |
| `EvolveECartesianECT` | E 更新乘 `1/ε̃`（Ex/Ez 用法向分量、Ey 用切向分量）；σ 项改 CN `alpha/beta` |
| B 侧 | 不动（B 是基本量）；`EvolveECTRho` 不动 |
| 稳定化 | **沿用 M2a 的 BCK 下限 + DM 冻结**；介质不引入新的小胞元机制（ε 大只收紧 dt，σ 已 CN） |
| 同源 | ε̃ 由同一 level set 构造（与 Lx/Lz/A_y 同源），避免 M2a 踩过的"不同源"坑 |

### 4.6 与 PEC 共存、PML、能量

- PEC+介质混合：PEC 侧 ε→∞ 不参与平均（照 M2a 掩膜跳过）；介质-PEC 界面用 ε_max 截断或解析取极限；
- PML：介质不进入 PML（沿用 M2a 的断言）；PML 内保持 Yee；
- 能量：无损 + 对称张量 ε̃ 时离散能量为正定（Kottke 的微扰理论保证到二阶）——加"能量漂移"验收行。

### 4.7 "嵌入 Riemann"的定位（重要澄清）

- FDTD 的界面条件（切向 E/H、法向 D/B 连续）在 Yee 网格上由 **§4.2 的平均规则自动满足**（这正是 Hwang–Cangellaris
  "有效 ε" 的内容）；**不需要也不应该**引入 Riemann 通量（FDTD 无通量）。
- "嵌入 Riemann"属于 **FV/DG 侧**：流体（M2a-2 的 PEC 镜像 Riemann，Colella 2006 同族）与 Gkeyll 风格 Maxwell-DG
  （architecture 决策 4 的 "COPY_B 思路"）。若将来把 Maxwell 也 FV 化（DG 二期），再谈。
- 文档建议：把 architecture 决策 4 的措辞更新为"场侧用共形等效 ε/μ（subpixel），流体侧用嵌入 Riemann"。

---

## 5. Fresnel R/T 验证阶梯（M2b 验收）

| 步 | 算例 | 验收 |
|---|---|---|
| V1 | 垂直入射介质平板（ε 扫描 1.5/2.25/4/16；TE 与 TM） | R/T 与 Fresnel 解析：\|R\| ≤1%、相位 ≤2°；分辨率收敛 ≥1.8 |
| V2 | 斜入射 15/30/45/60°，s/p 偏振 | 同上；Brewster 角处 R_p → 0（对无损耗介质） |
| V3 | 临界角/全反射（ε 内→外） | 全反射 \|R\|=1、相位符合；衰减场指数正确 |
| V4 | σ 扫描（有损耗平板，TE 垂直） | R/T 与复 Fresnel 解析；CN 稳定性（大 σ 不炸） |
| V5 | μ≠1 平板 | R/T 解析（Z=√(μ/ε) 阻抗） |
| V6 | 曲面介质（圆柱/透镜）+ 梯度标量场 ε(x,z) | 与 **Meep 对照**（同网格/同 level set） |
| V7 | PEC+介质混合（介质包 PEC） | R/T + 无热斑三行验收（全网格 max\|E\|、PEC 内、环带） |
| V8（二期） | 色散平板（Debye/Lorentz） | 透射谱 vs 解析/ gprMax 对照 |

判据与 M2a 同规：几何自检（规则域逐位退化）、能量漂移行、晚时监控行。

---

## 6. 实施顺序与工作量

| 步 | 内容 | 验收 | 工时 |
|---|---|---|---|
| M2b-1 | 分区常数 ε/μ/σ + CN-σ + 垂直入射平板 | V1（R/T ≤1%） | 1–1.5 天 |
| M2b-2 | 标量场 ε(x,z)/μ/σ + 等效系数（切向算术/法向调和，先标量版） | V2–V5 | 1 天 |
| M2b-3 | Kottke 各向异性张量 ε̃（level set 梯度 n̂ + 逐 cell 积分）+ 斜/曲面界面 | V2/V6/V7；与 Meep 对照 | 1.5–2 天 |
| M2b-4 | 色散 ADE（Debye/Drude/Lorentz） | V8 | 1–1.5 天 |
| M2b-5 | 文档定稿 + 基线重记 | 全套 | 0.5 天 |

顺序依赖：M2b-1 → M2b-2 → M2b-3 → M2b-4；M2b-1 尽量复用 `MacroscopicProperties` 与 `T_MacroAlgo::alpha/beta`。

---

## 7. 风险与开放问题

1. **ε̃ 与 AMReX 几何同源**：M2a 的教训（不同源 → O(h) 漂移）在此同样适用——ε̃ 必须由与 `Lx/Lz/A_y` 相同的
   level-set 积分构造（§4.5 同源行）。
2. **σ 的界面平均**：Kottke 的 σ̃ 规则在 cut cell 上需复算验证（V4 的衰减率对标会暴露问题）。
3. **PEC+介质界面的极限**：ε→∞ 的调和平均退化为 0/∞，需显式处理（掩膜跳过 + 极限截断）。
4. **大 ε 与 dt**：ε 大 → 波速低，全局 dt 不变时精度/稳定余量收紧；需要 ε_max 扫描验收。
5. **色散 + subpixel**：Meep 明说 subpixel 不适用色散介质（只取 ε∞）——M2b-4 的界面精度预期要降级，
   文档里写清"色散界面一阶"。
6. **能量/守恒**：无损 + 对称 ε̃ 时正定；有损/色散时按 gprMax 稳定性报告的做法记录能量趋势即可。

---

## 8. 参考文献（链接见 §2）

1. Yu & Mittra 2001（共形介质曲面）；2. Hwang & Cangellaris 2001（有效 ε，二阶）；
3. Kottke et al. 2008（各向异性微扰理论，核心公式）；4. Farjadpour et al. 2006（subpixel）；
5. Deinega & Valuev 2007（斜界面）；6. Taflove & Hagness 2005（ADE/CN-σ）；
7. Meep `Subpixel_Smoothing.md` 与 `src/anisotropic_averaging.cpp`（本机）；
8. WarpX `MacroscopicProperties`、`MacroscopicEvolveE.cpp`（本机）；
9. gprMax 材料稳定性报告（本机）。

## 附录：代码位置

| 件 | 位置 |
|---|---|
| Meep 有效张量公式 | `<local-ref>/meep/src/anisotropic_averaging.cpp:201`（公式）、`:100`（法向） |
| Meep 局限说明 | `<local-ref>/meep/doc/docs/Subpixel_Smoothing.md`（5 条限制） |
| WarpX 宏观属性 | `Source/FieldSolver/FiniteDifferenceSolver/MacroscopicProperties/` |
| WarpX σ CN 系数 | `Source/FieldSolver/FiniteDifferenceSolver/MacroscopicEvolveE.cpp:186-187`（`T_MacroAlgo::alpha/beta`） |
| ECT 内核接缝 | `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp`（`EvolveECartesianECT`）、`EmbeddedBoundaryInit.cpp`（`ComputeDualGeometry`） |

*本文是 M2b 的方法调研；M2a 场侧（G1–G4）与流体侧（M2a-2）的进度见 `C_SIDE_LS_ECT_PLAN.md` 与 `HANDOFF.md`。*
</content>
