# Cut-cell 方法文献调研报告（面向 M2a / M2b）

> 被调研文档：`<local-doc>/cutcell-research.md`
> （v2，130 行，mtime 2026-10-31 10:04；含 §7 决策与 §8 M2a 设计规范）
>
> 调研范围：共形 FDTD（PEC）、Enlarged Cell Technique、介质界面离散、EB 有限体积小胞元稳定化、AMR 回流、
> 以及用于验证的解析基准。文献全部经 Crossref 核对 DOI（2026-10-31）。
>
> 内部证据：`benchmarks/cases/c02_0_pec/ROOTCAUSE_ECT_2D.md`（2D XZ ECT 仅 TM）与本机代码（WarpX_fluid）。
>
> 结论速览：**文档的三项主要决策（E 侧 Dey–Mittra + ECT、B 侧复用 WarpX ECT、流体用 AMReX SRD）
> 都有明确的文献支持；但"E 侧与 B 侧对称"必须先做对偶几何（dual-grid）推导，不能直接复用
> `face_areas`；M2a 的 PEC 验证在 2D 只能做 TM，s 偏振需 3D 或补 2D TE。**

---

## 0. 决策 × 文献支持度

| 文档决策 | 支持度 | 关键文献 | 必须补的细节 |
|---|---|---|---|
| E 场：Dey–Mittra 共形 Ampere + ECT 面积扩展 | **强** | Dey–Mittra 1997；Xiao–Liu 2008；Jurgens–Taflove 1992 | **对偶面/对偶边**的几何与 `S_stab`；与 B 侧的 `Venl` 式通量累计配对 |
| B 场：复用 WarpX `EvolveBCartesianECT` | **强**（但有已知缺陷） | Xiao–Liu 2008；Benkler 2006；WarpX 代码 | 2D XZ 仅 TM（内部证据）；BCK 面问题（上游 issue #7206）|
| 流体：AMReX `EB_StateRedistribute` | **强** | Berger–Giuliani 2021；Giuliani 2022；Berger–Giuliani 2024 | 源项（Lorentz 力）的处理与守恒检验 |
| 介质 ε/μ/σ 推迟到 M2b | **合理** | Wang–Yin–Liu 2009；Farjadpour 2006；Kottke 2008；Deinega 2007 | M2b 三条路线（共形介质 / subpixel / 体积分数+Riemann）的取舍标准 |
| M2a 只做 PEC 闭环 | **合理** | Dey–Mittra 1997；Xiao–Liu 2008；Jurgens–Taflove 1992 | 验证需补：收敛阶、能量守恒、PEC 球/柱 Mie 基准、TE 偏振（3D）|

---

## 1. 调研方法

- **检索**：Crossref `query.bibliographic` 逐条核题名/卷期/页码/DOI；对关键文献做前向引用（谁引用了它）
  与后向引用（它引用了谁），形成三条方法谱系（§2）。
- **核对**：所有 DOI 均解析成功；卷期页码以 Crossref 记录为准（见 §8）。
- **范围**：只收与"PEC/介质 + cut-cell + 显式时域 + 小胞元"直接相关的文献；不展开到频域 FEM/DG。
- **内部证据**：WarpX 代码行号、2D TE 复现实验（见 ROOTCAUSE 报告），仅作"实现现实"标注，不冒充文献。

---

## 2. 方法谱系

### 2.1 共形 FDTD（PEC）：Jurgens–Taflove → Dey–Mittra → ECT → BCK

| 文献 | 贡献 | 对 M2a 的意义 |
|---|---|---|
| **Jurgens & Taflove 1992**（IEEE TAP 40:357–366） | 用"轮廓路径"（contour path）把弯曲 PEC 表面嵌入 Yee 网格，系统化共形 FDTD 的起点 | 说明"局部修改更新系数"可行；DM 在其基础上进一步改进 |
| **Dey & Mittra 1997**（IEEE MGWL 7:273–275） | 局部共形 FDTD（"DM 格式"）：cut 边用实际边长、cut 面用实际面积做 Faraday/Ampere 积分离散 | **M2a 的 E 侧模板**；给出逐边/逐面的系数修改方式 |
| **Dey & Mittra 1998**（MOTL 17:349） | 改进版（标题即 "modified"），针对原始格式的精度/稳定性 | 后续稳定化讨论（ECT/BCK）的起点 |
| **Xiao & Liu 2004**（IEEE MWCL 14:551–553） | "enlarged cells"：把小 cut 面从邻居借面积，避免时间步缩减 | ECT 的起点；WarpX 的 `area_mod`/`borrowing` 即此思想 |
| **Xiao & Liu 2005**（IEEE AP-S） | 把 enlarged cell 技术用于共形 FDTD 建模 PEC | 会议版，给出二维/三维的初步构造 |
| **Xiao & Liu 2008**（IEEE TAP 56:765–773） | **3-D ECT** 的完整构造 + 稳定性判据；对 E/H 两侧给出统一的"扩大面 + 借用"处理 | **M2a 的核心参考**：E 侧"与 B 侧对称"必须以该文的 dual 构造为准 |
| **Benkler, Chavannes & Kuster 2006**（IEEE TAP 54:1843–1849） | 另一种共形 PEC 格式 + 稳定性判据（"BCK"），用面积缩放稳定小面 | WarpX 用它做兜底（`ApplyBCKCorrection`）；E 侧也要有同样的 fallback |

**谱系结论**：DM 提供"局部共形"的离散模板；ECT 提供"小胞元稳定化"的机制；BCK 是另一种稳定化兜底。
三者不是互斥选项，而是**同一方案的不同层次**——DM 是离散，ECT/BCK 是稳定化。文档把"DM Ampere + ECT 扩展"
写在一起是正确的组合，但必须按 Xiao–Liu 2008 的**成对（E 与 H 两侧）**方式实现。

### 2.2 EB 有限体积的小胞元稳定化：h-box → 通量重分配 → 单元合并 → SRD

| 文献 | 方法 | 特点 | 对 M2a 流体的意义 |
|---|---|---|---|
| **Berger 2003**（SIAM J. Numer. Anal. 41:893–918）；**Helzel, Berger & LeVeque 2005**（SISC 26:785–809） | **h-box**：把 cut cell 与邻居合并成"h-box"，在 h-box 上做一维/旋转 Riemann 求解 | 可用全胞元 CFL；二阶；实现复杂、需要构造 h-box | 备选（若 SRD 在强激波下不满意） |
| **Pember et al. 1995**（JCP 120:278–304）；**Colella et al. 2006**（JCP 211:347–366） | **flux redistribution**：小胞元的超额通量按权重分给邻居 | 稳定、实现相对简单；边界附近精度/阶数会受影响 | 备选；AMReX 的 `apply_flux_redistribution` 即此 |
| **Berger 2015**（Appl. Numer. Math. 96:180–186） | **cell merging** 的稳定性分析 | 稳定，但改变网格拓扑、实现重 | 备选 |
| **Berger & Giuliani 2021**（JCP 428:109820） | **state redistribution (SRD)**：每步后把 cut cell 的"差量"按体积分数重分配给邻居 | 守恒、保线性、**不降阶**、显式、按全胞元体积取 dt | **文档选它的理由充分** |
| **Giuliani et al. 2022**（JCP 464:111305） | **weighted SRD**：加权最小二乘改善小胞元附近精度 | 在 SRD 基础上更稳、更准 | 若基础 SRD 在曲面附近精度不够，升级到 weighted |
| **Berger & Giuliani 2024**（SISC 46(5):A2848–A2873） | 可证稳定的 weighted SRD（单调、TVD、GKS 稳定） | 理论保证最强 | 实现时可作为正确性依据 |

**结论**：文档选 SRD 与文献主流一致，且 AMReX 已有实现（`AMReX_EB_StateRedistribute` /
`AMReX_EB_Redistribution`）。h-box 与 flux redistribution 作为对照/回退保留在阅读清单里即可。

### 2.3 介质界面：subpixel vs 共形介质

| 文献 | 方法 | 特点 |
|---|---|---|
| **Farjadpour et al. 2006**（Opt. Lett. 31:2972） | Meep 的 subpixel smoothing：法向调和、切向算术，构造各向异性 ε⁻¹ | 论文报告对曲面界面接近二阶、实现简单；M2b 的首选基线 |
| **Kottke, Farjadpour & Johnson 2008**（PRE 77:036611） | 各向异性介质界面的微扰理论，解释 subpixel 为什么有效、误差从哪来 | 给出 subpixel 的理论依据与改进方向 |
| **Deinega & Valuev 2007**（Opt. Lett. 32:3429） | 把 subpixel 推广到有耗/色散介质 | 对应文档里的 σ（有耗）路线 |
| **Wang, Yin & Liu 2009**（Electron. Lett. 45:146–147） | (2,4) 高阶共形介质界面技术 | 高阶、复杂；可作为 M2b 的"高精度对照" |

**结论**：M2b 建议先做 subpixel（简单、二阶、文献充分），再与共形介质比较；不要一上来就做高阶共形。

### 2.4 替代路线（用于回答"为什么不选"）

| 文献 | 路线 | 不选的理由 |
|---|---|---|
| Shankar, Mohammadian & Hall 1990（Electromagnetics 10:127） | FVTD（有限体积时域）求解 Maxwell | 需要非结构/切割网格与 Riemann 通量，和 WarpX 的 Yee 场布局不一致；迁移成本高 |
| Ingram et al. 2003（Math. Comput. Simul. 61:561） | Cartesian cut-cell 方法综述（CFD 视角） | 提供小胞元问题的全景；作为背景阅读 |
| Popescu, Vedder & Shyy 2008（IJNMF 56:1787–1818）* | 高阶 Cartesian cut-cell 波传播（OPC-fv） | 高阶、适合声学；与 FDTD 交错网格不直接兼容 |

*该条 DOI 10.1002/fld.1517，用于说明"高阶路线"存在，但不在 M2a 范围内。

---

## 3. 逐决策评估

### D1　E 场：Dey–Mittra 共形 Ampere + ECT 面积扩展

**文献怎么说**
- DM 1997 的 Ampere 更新是"对偶面"上的积分：E 在边上，更新用绕该边的 B 环量除以**该边的对偶面面积**；
  Faraday 更新则用 cut 边的实际长度。二者是一对（curl 的对偶）。
- ECT（Xiao–Liu 2004/2008）把"面积太小"的面抬高到稳定阈值 `S_stab`，并从邻居借面积；借来的通量必须
  在邻居处**记账**（WarpX 用 `Venl` 累加实现），否则离散守恒/一致性被破坏。
- 2008 年的 3D ECT 对 E 与 H **两侧**都定义了扩展，并给出稳定性判据；这正是"与 B 场对称"的出处。

**文档的正确之处**：把"DM Ampere"与"ECT 扩展"组合，方向正确；"注意 E/B 面错位半格需重新取几何量"
这句话也点到了要害。

**文档缺的（必须在实现前补）**
1. **对偶几何**：E 边的对偶面面积**不是** WarpX 的 `face_areas`（后者是 B 面）。需要为每条 E 边构造
   对偶面（面积 = 四个相邻 B 面的四分之一面积之和，或按 cut 几何直接积分），并给出对偶边长度。
2. **对偶版的 `S_stab` 与 borrowing**：Xiao–Liu 2008 的判据要按对偶面重推，不能直接套用 WarpX
   现有的面版 `ComputeSStab`（`WarpXFaceExtensions.cpp:88`）。
3. **通量记账**：B 侧用 `Venl` 把借来的环量累加；E 侧需要一个对应的"借来的 B 通量"累加器，
   否则扩展面积上的更新与邻居的更新不闭合。
4. **BCK fallback**：B 侧在"借不到面积"时退到 BCK（`ApplyBCKCorrection`）；E 侧要有同样的分支。

**建议**：保留 D1，但把它拆成 **M2a-1a（对偶几何 + 判据）** 与 **M2a-1b（更新内核 + 记账）** 两个子任务，
先在一个 1D/2D 的斜 PEC 平面上把"对偶面面积 = 解析值"验证到机器精度，再进 3D。

### D2　B 场：复用 WarpX `EvolveBCartesianECT`

**文献怎么说**：WarpX 的 B 侧实现对应 Xiao–Liu 2008 的 H 更新（面积扩展 + borrowing + BCK 兜底），
三分支（extended / available+bck / intruded）与 `Venl` 记账在代码里都能对上（`EvolveB.cpp:285-395`）。

**文档没写、但会影响 M2a 的两件事**
1. **2D XZ 只实现了 TM**：`EvolveECTRho.cpp` 在 XZ 下不计算 Rhox/Rhoz，`EvolveB.cpp:261` 的循环
   只到 `AMREX_SPACEDIM`（2D 只更新 Bx、By），**Bz 永不更新、Bx 用 Rho=0 更新**。受控实验：
   2D ECT + TE（Ey+Bx）→ B 冻结、E 线性增长；2D ECT + TM（Ex+By）→ 正常；3D ECT + TE → 正常；
   2D Yee + TE → 正常。详见 `ROOTCAUSE_ECT_2D.md`。
   **含义**：M2a 的 PEC 验证在 2D 只能做 TM；s 偏振/ TE 模要么用 3D，要么先补 2D TE。
2. **上游 issue #7206**（"ECT: B is never updated on faces stabilized by the BCK correction"）
   提示 BCK 面在某个版本上存在更新缺失；虽然当前代码看起来 BCK 分支有 `B -= dt*Rho`，
   但在做 E 侧对称实现前应把这条 issue 的复现条件核对一遍，避免把 bug 复制到 E 侧。

**建议**：D2 保持；在 M2a 的"场侧"工作里加一条**前置检查**：2D/3D、TM/TE 的组合矩阵必须先跑通
（用 §5 的基准），再把 E 侧接上。

### D3　流体：AMReX SRD

**文献怎么说**：Berger–Giuliani 2021 明确 SRD 的三个卖点——守恒、保线性（不降阶）、按全胞元体积取 dt；
2022 的 weighted 版改善小胞元附近精度；2024 给出可证稳定（单调/TVD/GKS）。AMReX 同时提供
`single_level_redistribute`（基础）与 `single_level_weighted_redistribute`（加权），以及 `apply_flux_redistribution`。

**文档的正确之处**：选 SRD、复用 AMReX，都是文献支持的最省力路径。

**文档缺的**
1. **源项**：流体矩方程带 Lorentz 力、压力/热流源项；SRD 的文献针对守恒型双曲系统（无源项），
   带源项时的守恒/精度没有现成结论。建议**算子分裂**：先做 FV 双曲步 + SRD，再加源项，
   并单独验证"源项步不破坏守恒"。
2. **边界条件**：PEC 面上流体"反射（无穿透）"需要 ghost-cell/Riemann 处理（Colella 2006 的 EB
   边界做法），文档只写了"反射"，没写离散方式。
3. **与场耦合**：SRD 之后沉积的电荷/电流与场的守恒（∇·E−ρ）需要单独检查；文献没有覆盖
   "SRD + Maxwell 源项"的组合。

**建议**：D3 保持，但在 M2a-2 里把"源项 + 边界 + 场耦合"单列为验收项，不要只验"守恒误差"。

### D4　介质推迟到 M2b

**文献怎么说**：PEC 的 DM/ECT 是成熟方法；介质界面没有同等程度的"标准答案"：
- subpixel（Farjadpour 2006；Kottke 2008；Deinega 2007）：简单、二阶，但对"切割网格上的 FV 流体 + Maxwell"
  需要把 ε 的有效张量接到 Yee 更新里；
- 共形介质（Wang 2009）：高阶、复杂；
- 体积分数 + 界面 Riemann（Colella 2006 一类）：适合 FV，但要把 Maxwell 写成守恒形式并处理
  ε/μ 的界面条件。

**建议**：M2b 先做 subpixel（标量法向投影），用 Fresnel/Brewster 与介质平板基准验收；通过后再评估
是否值得上共形介质。文档的"三条路线"应加一条**选择顺序**：subpixel → 共形介质 → 体积分数+Riemann。

### D5　M2a 只做 PEC

**文献怎么说**：DM/ECT/BCK 都是 PEC 方法；PEC 的验证基准在文献里很成熟：PEC 球/柱的 RCS 与 Mie 级数、
矩形/圆柱腔模式、波导与绕射。Jurgens–Taflove 与 DM 都用 RCS 曲线做验证；Xiao–Liu 2008 也用 RCS/散射。

**文档的验证设计（§8.4）可以，但不够**
- 缺**收敛阶**：至少 3 个分辨率，报告误差随 Δx 的斜率；
- 缺**能量/长时间稳定性**：共形格式的一个卖点是长时间稳定，必须测 100+ 周期；
- 缺**解析基准**：PEC 球/柱的 RCS 与 Mie 级数对比（不是只看"反射精度"）；
- 缺**偏振矩阵**：2D 只能 TM（见 D2），s 偏振进 3D；
- 缺**空腔模式**：`Examples/Tests/embedded_boundary_rotated_cube` 已有 2D/3D 腔模式测试，可直接复用。

---

## 4. 实现细节的文献依据（可直接落地）

| 细节 | 文献依据 | 落地建议 |
|---|---|---|
| DM 的 E/H 更新系数 | Dey–Mittra 1997；Jurgens–Taflove 1992 | 先在规则网格上验证"退化回 Yee"（对偶面=Δ²、边=Δ 时逐位一致） |
| `S_stab` 的定义 | Xiao–Liu 2008；WarpX `ComputeSStab`（0.5·max{邻边×间距}） | E 侧要推对偶版；用"规则网格恒等于 Δ²"做单测 |
| 面积借用与记账 | Xiao–Liu 2008；WarpX `Venl`/`borrowing` | E 侧复制 `Venl` 模式；单测"借出的通量 = 借入的通量" |
| BCK 兜底 | Benkler et al. 2006；WarpX `ApplyBCKCorrection` | E 侧同样要有；核对上游 issue #7206 |
| SRD 步骤 | Berger–Giuliani 2021（算法）；2022（加权）；2024（证明） | 直接用 AMReX API；单测"均匀流场经 SRD 逐位不变" |
| 流体 EB 边界 | Colella et al. 2006；Pember et al. 1995 | 用 ghost-cell + Riemann；PEC 壁面反射用镜像 |
| AMR 回流 | Berger–Colella 1989；AMReX `EBFluxRegister` | 留到 M4；M2a 固定单层 |

---

## 5. 验证基准（文献用过的，可直接采用）

| 基准 | 解析/参考解 | 文献 |
|---|---|---|
| PEC 球 RCS | Mie 级数（PEC 极限） | Jurgens–Taflove 1992；Dey–Mittra 1997 |
| PEC 柱 RCS（2D） | 柱 Mie 级数 | DM/ECT 文献；WarpX 2D 测试（仅 TM） |
| 矩形/圆柱腔模式 | 解析本征频率 | WarpX `embedded_boundary_rotated_cube` |
| 斜 PEC 平面反射 | $r=-1$、相位反相 | 本项目 C2-0（注意 2D 偏振限制） |
| 介质平板 Fresnel / Brewster | 复 Fresnel | M2b 用；见 `EM_EB_MHD_VALIDATION_10CASES_v2.md` C1/C3 |
| 流体绕圆柱/激波撞壁 | 实验/高分辨率参考 | Pember 1995；Colella 2006；Helzel 2005 |

**建议的 M2a 验收（补充文档 §8.4）**
1. 3D 斜 PEC 平面（s 与 p 两种偏振）：$|r+1|\le1\%$、相位 $\le2^\circ$、$\epsilon_{\rm spur}\le10^{-3}$；
2. 3D PEC 球：RCS 与 Mie 级数比，$\le2\%$（$ka\ge1$）；
3. 2D 圆柱（TM）：与柱 Mie 比，$\le3\%$；
4. 腔模式：本征频率 $\le0.5\%$（直接复用 WarpX 测试）；
5. 收敛阶：$\ge1.5$（3 个分辨率）；
6. 能量：100 周期漂移 $\le10^{-3}$；
7. 流体：SRD 前后总质量/动量/能量守恒 $\le10^{-10}$（无源项步），带源项时 $\le10^{-6}$。

---

## 6. 风险与未决问题

| 风险 | 来源 | 缓解 |
|---|---|---|
| E 侧对偶几何推导错误 | 文献只给 3D 构造；2D 需自己推 | 先做"规则网格退化"单测；再做斜平面解析对偶面积对比 |
| 2D ECT 只支持 TM | 内部证据（ROOTCAUSE） | M2a 的 2D 验证限定 TM；TE 用 3D，或补 2D TE |
| BCK 面更新问题 | 上游 issue #7206 | 实现前复现并确认；E 侧不要复制错误 |
| SRD + 源项/场耦合无文献 | 文献只覆盖守恒双曲 + 无源 | 算子分裂 + 单独守恒检验；把结论写进报告而非假设 |
| 介质（M2b）无标准答案 | 三条路线各有取舍 | 先 subpixel，再比较；Fresnel 基准把关 |
| 高阶/低阶混淆 | 文档未规定阶数目标 | 明确 M2a 目标为二阶（DM/ECT 的典型精度），不追高阶 |

---

## 7. 阅读清单（按优先级）

**P0（实现前必读）**
1. Dey & Mittra 1997（DM 格式）
2. Xiao & Liu 2008（3-D ECT 完整构造与判据）
3. Berger & Giuliani 2021（SRD 算法）
4. Farjadpour et al. 2006（subpixel；M2b 用）
5. AMReX EB 文档 + `AMReX_EB_StateRedistribute` 源码

**P1（设计与验证）**
6. Jurgens & Taflove 1992（共形 FDTD 的起点）
7. Xiao & Liu 2004（enlarged cells 的动机）
8. Benkler et al. 2006（BCK 兜底）
9. Wang, Yin & Liu 2009（介质共形）
10. Colella et al. 2006；Pember et al. 1995（EB 边界与通量重分配）
11. Giuliani et al. 2022（weighted SRD）

**P2（对照与背景）**
12. Helzel, Berger & LeVeque 2005；Berger 2003（h-box）
13. Berger 2015（cell merging 稳定性）
14. Berger & Giuliani 2024（可证稳定 SRD）
15. Kottke et al. 2008；Deinega & Valuev 2007（subpixel 理论/有耗）
16. Ingram et al. 2003（cut-cell 综述）
17. Shankar et al. 1990（FVTD 替代路线）

---

## 8. 参考文献（DOI 均经 Crossref 核对，2026-10-31）

1. K. S. Jurgens, A. Taflove, *Finite-difference time-domain modeling of curved surfaces*,
   IEEE Trans. Antennas Propag. 40(4):357–366 (1992).
   DOI: [10.1109/8.138836](https://doi.org/10.1109/8.138836)
2. S. Dey, R. Mittra, *A locally conformal finite-difference time-domain (FDTD) algorithm for
   modeling three-dimensional perfectly conducting objects*, IEEE Microw. Guided Wave Lett.
   7(9):273–275 (1997). DOI: [10.1109/75.622536](https://doi.org/10.1109/75.622536)
3. S. Dey, R. Mittra, *A modified locally conformal finite-difference time-domain algorithm for
   modeling three-dimensional perfectly conducting objects*, Microw. Opt. Technol. Lett.
   17(6):349–352 (1998). DOI: [10.1002/(SICI)1098-2760(19980420)17:6<349::AID-MOP4>3.0.CO;2-H](https://doi.org/10.1002/%28SICI%291098-2760%2819980420%2917:6%3C349::AID-MOP4%3E3.0.CO;2-H)
4. T. Xiao, Q. H. Liu, *Enlarged cells for the conformal FDTD method to avoid the time step
   reduction*, IEEE Microw. Wirel. Components Lett. 14(12):551–553 (2004).
   DOI: [10.1109/LMWC.2004.837384](https://doi.org/10.1109/LMWC.2004.837384)
5. T. Xiao, Q. H. Liu, *An enlarged cell technique for the conformal FDTD method to model
   perfectly conducting objects*, IEEE AP-S Int. Symp. (2005).
   DOI: [10.1109/APS.2005.1551259](https://doi.org/10.1109/APS.2005.1551259)
6. T. Xiao, Q. H. Liu, *A 3-D enlarged cell technique (ECT) for the conformal FDTD method*,
   IEEE Trans. Antennas Propag. 56(3):765–773 (2008).
   DOI: [10.1109/TAP.2008.916876](https://doi.org/10.1109/TAP.2008.916876)
7. S. Benkler, N. Chavannes, N. Kuster, *A new 3-D conformal PEC FDTD scheme with user-defined
   geometric precision and derived stability criterion*, IEEE Trans. Antennas Propag.
   54(6):1843–1849 (2006). DOI: [10.1109/TAP.2006.875909](https://doi.org/10.1109/TAP.2006.875909)
8. J. Wang, W.-Y. Yin, Q.-H. Liu, *FDTD (2, 4)-compatible conformal technique for treatment of
   dielectric surfaces*, Electron. Lett. 45(3):146–147 (2009).
   DOI: [10.1049/el:20093215](https://doi.org/10.1049/el:20093215)
9. A. Farjadpour *et al.*, *Improving accuracy by subpixel smoothing in the finite-difference
   time domain*, Opt. Lett. 31(20):2972–2974 (2006).
   DOI: [10.1364/OL.31.002972](https://doi.org/10.1364/OL.31.002972)
10. A. Kottke, A. Farjadpour, S. G. Johnson, *Perturbation theory for anisotropic dielectric
    interfaces, and application to subpixel smoothing of discretized numerical methods*,
    Phys. Rev. E 77:036611 (2008). DOI: [10.1103/PhysRevE.77.036611](https://doi.org/10.1103/PhysRevE.77.036611)
11. A. Deinega, I. Valuev, *Subpixel smoothing for conductive and dispersive media in the
    finite-difference time-domain method*, Opt. Lett. 32(23):3429–3431 (2007).
    DOI: [10.1364/OL.32.003429](https://doi.org/10.1364/OL.32.003429)
12. C. Helzel, M. J. Berger, R. J. LeVeque, *A high-resolution rotated grid method for
    conservation laws with embedded geometries*, SIAM J. Sci. Comput. 26(3):785–809 (2005).
    DOI: [10.1137/S106482750343028X](https://doi.org/10.1137/S106482750343028X)
13. M. J. Berger, *H-box methods for the approximation of hyperbolic conservation laws on
    irregular grids*, SIAM J. Numer. Anal. 41:893–918 (2003).
    DOI: [10.1137/S0036142902405394](https://doi.org/10.1137/S0036142902405394)
14. R. B. Pember, J. B. Bell, P. Colella, W. Y. Crutchfield, M. L. Welcome, *An adaptive
    Cartesian grid method for unsteady compressible flow in irregular regions*,
    J. Comput. Phys. 120(2):278–304 (1995).
    DOI: [10.1006/jcph.1995.1165](https://doi.org/10.1006/jcph.1995.1165)
15. P. Colella, D. T. Graves, B. J. Keen, D. Modiano, *A Cartesian grid embedded boundary method
    for hyperbolic conservation laws*, J. Comput. Phys. 211(1):347–366 (2006).
    DOI: [10.1016/j.jcp.2005.05.026](https://doi.org/10.1016/j.jcp.2005.05.026)
16. M. J. Berger, *A note on the stability of cut cells and cell merging*,
    Appl. Numer. Math. 96:180–186 (2015).
    DOI: [10.1016/j.apnum.2015.05.003](https://doi.org/10.1016/j.apnum.2015.05.003)
17. M. J. Berger, A. Giuliani, *A state redistribution algorithm for finite volume schemes on
    cut cell meshes*, J. Comput. Phys. 428:109820 (2021).
    DOI: [10.1016/j.jcp.2020.109820](https://doi.org/10.1016/j.jcp.2020.109820)
18. A. Giuliani, A. S. Almgren, J. B. Bell, M. J. Berger, M. T. Henry de Frahan, D. Rangarajan,
    *A weighted state redistribution algorithm for embedded boundary grids*,
    J. Comput. Phys. 464:111305 (2022).
    DOI: [10.1016/j.jcp.2022.111305](https://doi.org/10.1016/j.jcp.2022.111305)
19. M. J. Berger, A. Giuliani, *A new provably stable weighted state redistribution algorithm*,
    SIAM J. Sci. Comput. 46(5):A2848–A2873 (2024).
    DOI: [10.1137/23M1597484](https://doi.org/10.1137/23M1597484)
20. D. M. Ingram, D. M. Causon, C. G. Mingham, *Developments in Cartesian cut cell methods*,
    Math. Comput. Simul. 61:561–572 (2003).
    DOI: [10.1016/S0378-4754(02)00107-6](https://doi.org/10.1016/S0378-4754(02)00107-6)
21. V. Shankar, A. H. Mohammadian, W. F. Hall, *A time-domain, finite-volume treatment for the
    Maxwell equations*, Electromagnetics 10(1–2):127–145 (1990).
    DOI: [10.1080/02726349008908232](https://doi.org/10.1080/02726349008908232)
22. W. Zhang *et al.*, *AMReX: a framework for block-structured adaptive mesh refinement*,
    J. Open Source Softw. 4(37):1370 (2019).
    DOI: [10.21105/joss.01370](https://doi.org/10.21105/joss.01370)
23. W. Zhang *et al.*, *AMReX: Block-structured adaptive mesh refinement for multiphysics
    applications*, Int. J. High Perform. Comput. Appl. 35(6):512–526 (2021).
    DOI: [10.1177/10943420211022811](https://doi.org/10.1177/10943420211022811)
24. M. J. Berger, J. Oliger, *Adaptive mesh refinement for hyperbolic partial differential
    equations*, J. Comput. Phys. 53(3):484–512 (1984).
    DOI: [10.1016/0021-9991(84)90073-1](https://doi.org/10.1016/0021-9991(84)90073-1)
25. M. J. Berger, P. Colella, *Local adaptive mesh refinement for shock hydrodynamics*,
    J. Comput. Phys. 82(1):64–84 (1989).
    DOI: [10.1016/0021-9991(89)90035-1](https://doi.org/10.1016/0021-9991(89)90035-1)
26. M. J. Berger, *On conservation at grid interfaces*, SIAM J. Numer. Anal. 24(5):967–984 (1987).
    DOI: [10.1137/0724063](https://doi.org/10.1137/0724063)

30. M. Popescu, R. Vedder, W. Shyy, *A finite volume-based high-order, Cartesian cut-cell
    method for wave propagation*, Int. J. Numer. Methods Fluids 56(10):1787–1818 (2008).
    DOI: [10.1002/fld.1517](https://doi.org/10.1002/fld.1517)

**教科书**
27. A. Taflove, S. C. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain
    Method*, 3rd ed., Artech House (2005)（FDTD 与共形/稳定性标准参考）。
28. R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge Univ. Press (2002)。
29. C. K. Birdsall, A. B. Langdon, *Plasma Physics via Computer Simulation*（PIC 噪声/色散，作背景）。

---

## 附录 A　文档论点 → 文献对照

| `cutcell-research.md` 的论点 | 文献判定 | 依据 |
|---|---|---|
| ECT = enlarged cell technique，避免时间步缩减 | ✅ | Xiao–Liu 2004/2008 |
| BCK 作为兜底 | ✅ | Benkler 2006；WarpX `ApplyBCKCorrection` |
| E 场仅掩膜、需共形 Ampere | ✅ | Dey–Mittra 1997；Jurgens–Taflove 1992 |
| 小胞元问题可用 SRD/通量重分配/h-box 解决 | ✅ | Berger–Giuliani 2021；Colella 2006；Helzel 2005 |
| 介质：subpixel / 共形介质 / 体积分数+Riemann | ✅（三路线都存在） | Farjadpour 2006；Wang 2009；Colella 2006 |
| "与 B 场对称"即可复用 ECT | ⚠️ 需对偶几何 | Xiao–Liu 2008 的 dual 构造；本项目 ROOTCAUSE 报告 |
| 2D 也能做 s 偏振 | ❌ 当前实现不支持 | 内部证据（2D XZ ECT 仅 TM） |

## 附录 B　检索记录（可复现）

```text
Crossref: query.bibliographic="locally conformal FDTD"            -> 10.1109/75.622536
Crossref: query.bibliographic="enlarged cell technique Xiao Liu"  -> 10.1109/TAP.2008.916876
Crossref: query.bibliographic="conformal PEC FDTD Benkler"        -> 10.1109/TAP.2006.875909
Crossref: query.bibliographic="state redistribution cut cell"     -> 10.1016/j.jcp.2020.109820
Crossref: query.bibliographic="weighted state redistribution"     -> 10.1016/j.jcp.2022.111305
Crossref: query.bibliographic="h-box irregular grids Berger"      -> 10.1137/S0036142902405394
Crossref: query.bibliographic="Cartesian grid embedded boundary"  -> 10.1016/j.jcp.2005.05.026
Crossref: query.bibliographic="subpixel smoothing Meep"           -> 10.1364/OL.31.002972
Crossref: query.bibliographic="anisotropic dielectric perturbation"-> 10.1103/PhysRevE.77.036611
```

> 本报告只做文献调研与决策评估，不含任何仿真结果；内部实验数据见 `ROOTCAUSE_ECT_2D.md`。
</content>
