# C 侧调研（第二版）：DESIGN2 的物理/代码核对、新发现、疑问清单与开工条件

> 日期：2026-10-31 · 状态：**调研更新，待实现**
> 对象：`C_SIDE_DESIGN2.md`（345 行，替代第一版 §3–§5）
> 定位：本文之于 `C_SIDE_DESIGN2.md`，等价于 `C_SIDE_RESEARCH.md` 之于第一版设计。
> 内容分四类：**问题（已确认）**、**疑问/存疑/质疑（带依据与置信度）**、**建议（带理由与替代方案）**、
> **文献/教科书/开源/社区调研**。第一版调研仍成立的文献主线不重复，只写增量。
>
> 结论速览：**DESIGN2 已修正第一版全部 4 个 P0，三条内核公式与源码逐位退化成立，可以开工；
> 但开工前必须先处理 2 个 P1（方案 A 半格偏移、PML box 未排除）与若干 P2；另有 16 条疑问/存疑
> 需要边做边验证（§5），其中 3 条（米氏能否被 BCK-only 根治、共形 E+PML 稳定性、TM 基线变化方向）
> 直接决定里程碑成败。**
> 底层物理：对偶配对的能量守恒要求两侧几何"同源"；混合取法只能保证一阶配对（§3）。
> 开源证据：Meep 的 subpixel 是 level-set 解析填充（与方案 C 同构）且**不适用于 PEC/有损材料**；
> gprMax 的 SIBC 稳定性报告把"小胞元 cut 面积可无限缩小"列为共形 FDTD 的公认问题，并给出
> "电对偶面积 ¼/½/¾"的同款构造先例（§4.3）。

---

## 0. 摘要

| 维度 | 结论 |
|---|---|
| DESIGN2 的 4 个 P0 修正 | ✅ 全部落实：`Ex`/`Ez` 分母 = `Lz`/`Lx`；`Ey` 单边中点值、无平均；退化算术自洽；存储拆成三套 |
| 三条内核公式 | ✅ 与 `EvolveE.cpp:167-199` 在规则域**逐位**一致（§1 复算） |
| B 侧对应表 | ✅ 与 fork 的 `ComputeFaceAreas` 一致；AMReX `areafrac` 的 lo 面约定已从源码确认 |
| **方案 A 半格偏移** | ❌ P1：`Lx`/`Lz` 公式与自身定义差半格（§2.1），需修正或改用方案 B |
| **PML box 未排除** | ❌ P1：`EvolveECartesian` 跑全网格、`EvolveEPML` 是叠加修正（§2.2） |
| 自检 F-2 空洞 / E-1 与稳定化冲突 | ⚠️ P2（§2.6/§2.7） |
| 疑问/存疑清单 | 16 条，含 3 条里程碑级风险（§5） |
| 文献/教科书/开源/社区 | 新增 Railton–Schneider、Gedney、Sullivan、Kuo–Kuo、Esirkepov；Meep/gprMax 源码级证据；WarpX #7206、Meep #1064（§4） |

---

## 1. DESIGN2 逐条核对（源码 + AMReX + 推导）

| # | DESIGN2 主张 | 核对 | 判定 |
|---|---|---|---|
| 1 | `Bx(i,j+½)` 在 cell 的 x-lo 面中心；`By` 在 cell 中心；`Bz(i+½,j)` 在 z-lo 面中心（§1.1） | 2D 交错：`Bx(1,0)`、`By(0,0)`、`Bz(0,1)`（`WarpX.cpp:2459-2474`） | ✅ |
| 2 | `face_areas[0]←getAreaFrac()[0]`、`[1]←getVolFrac()`、`[2]←getAreaFrac()[1]`（§1.1） | fork `ComputeFaceAreas` XZ 分支（`EmbeddedBoundaryInit.cpp:609-659`）正是如此 | ✅ |
| 3 | AMReX `getAreaFrac()[idim]` 是 **lo 侧**面（§1.5 陷阱注） | AMReX 源码 `AMReX_EBDataCollection.cpp:60-64`：`faceba = convert(a_ba, IntVect::TheDimensionVector(idim))` → nodal idim = 面的 lo 侧索引 | ✅（lo 面） |
| 4 | `A_y`/`Lx`/`Lz` 的位置与 `IntVect`（§1.1 表） | `A_y` 与 `Ey`(1,1) 同位；`Lx` 与 `Bx`(1,0) 同位；`Lz` 与 `Bz`(0,1) 同位 | ✅；注意 2D 构建里 `IntVect` 只有 2 分量，D-10 的"第三分量显式给 0"只在 3D 构建有含义 |
| 5 | 三条内核公式（§1.2） | 代入 `Lx=dx、Lz=dz、A_y=dxdz` 后与 `EvolveE.cpp:167-199` 逐项相等（含符号） | ✅ |
| 6 | `A_y = ¼ dx dz Σ v`（§1.4） | 四象限 ¼ 只属于节点量 `Ey` | ✅ |
| 7 | **方案 A**（§1.5） | **与 §1.5 自己的定义差半格**，见 §2.1 | ❌ P1 |
| 8 | 方案 B（§1.5） | 窗口正确：线段正好穿过这两个 cell 的中心线 | ✅（一阶，近竖直/近水平截面精确） |
| 9 | 方案 C（level set 采样）成本 <1 s（§1.5） | 量级合理；但 φ 符号要按项目约定（正值=固体），DESIGN2 写的 `[φ>0]=开放` 方向反了；"<1 s"是估计不是承诺（§5 Q13） | ⚠️ P2 |
| 10 | 开关三层退化（§2.1） | C-0/C-1a/C-1b 合理 | ✅；伪代码枚举名应为 `ElectromagneticSolverAlgo::ECT`，`WARPX_DIM_XZ` 是编译期宏、需 `#ifdef` 包裹 |
| 11 | covered 点判据 `<= 0` 并写 0（§2.3.2） | AMReX 对 covered 给精确 0 | ✅ |
| 12 | 掩膜新语义（§2.3.3） | 方向正确；但 `ly==0`（任一相邻边覆盖→冻结，`EmbeddedBoundaryInit.cpp:564-590`+`:258-292`）比 `A_y<=0` 更严 → 会继续冻结"cut 但开放"的节点，必须同步改 | ⚠️ P2（§2.4） |
| 13 | PML：新内核只作用 non-PML 区（§2.3.4） | **机制表述错**（§2.2） | ❌ P1 |
| 14 | J 项：非零即 abort（§2.3.5） | 比静默置零安全；建议同时检查"粒子数/电流源是否激活" | ✅（可接受） |
| 15 | 内存/生命周期/checkpoint（§2.4） | 三套 MultiFab nGrow=1；`Fields.H:133-135` 与 `:155-156` 需加三项 | ✅ |
| 16 | 稳定化：只抬分母副本、分子真实（§2.5） | 与 B 侧 BCK 对称（`Rho = 环量/S_enlarged`）；但 `Lx` 双重角色需要两套读数，且抬分母破坏对偶配对 → E-1 目标分档 | ✅ 方向对；⚠️ §3.3 |
| 17 | 自检 G-1/G-2/G-3/F-1/F-2/E-1（§3） | G-1/G-2/G-3/F-1 ✅；**F-2 空洞**（§2.6）；**E-1 冲突**（§2.7） | ⚠️ |
| 18 | S0'–S11'（§4） | 顺序合理；S1a' 应同时用解析平面验证"方案 B vs 修正后方案 A" | ✅（微调 §8） |
| 19 | 3D 推广表（§5） | 结构正确（对偶面四象限）；分子权重标签沿用易错命名（`Ex` 应为 `By×Ly、Bz×Lz`）；3D 同样有 lo/hi 半格问题 | ⚠️ P3 |
| 20 | 工作量 810 行 / 7–9 天（§6） | 合理；加"方案 A 修复 + PML 排除"约 0.5 天 | ✅ |

---

## 2. 本次新发现（开工前必读）

### 2.1【P1】方案 A 与它自己的定义差半格

DESIGN2 §1.5 的定义：

```
Lx(i,j+½) = 格心线 z=z_{j+½} 上、x∈[x_i−dx/2, x_i+dx/2] 的开放长度   // 跨 cell(i−1,j) 与 cell(i,j)
Lz(i+½,j) = 格心线 x=x_{i+½} 上、z∈[z_j−dz/2, z_j+dz/2] 的开放长度   // 跨 cell(i,j−1) 与 cell(i,j)
```

方案 A 的公式：

```
Lx(i,j+½) = ½ dx [a_z(i,j) + a_z(i,j+1)]    // 给出 cell(i,j) 的中心线 x∈[x_i,x_i+dx]
Lz(i+½,j) = ½ dz [a_x(i,j) + a_x(i+1,j)]    // 给出 cell(i,j) 的中心线 z∈[z_j,z_j+dz]
```

- `Lx` 公式给出的是 **cell (i,j) 的水平中心线**（`x∈[x_i,x_i+dx]`），定义要的是 `x∈[x_i−dx/2,x_i+dx/2]`——**x 向偏 dx/2**。
- `Lz` 同理在 **z 向偏 dz/2**。
- 反例：竖直 PEC 平面恰过 `x_i` 时，正确开放长度 = `dx/2`，公式给 `dx`。
- 根源：`Bx` 在 x 向是 nodal（对偶段跨两个 cell），`Bz` 在 z 向是 nodal；方案 A 只取了单个 cell 的两个面。

**修正（四面平均）**：

```
Lx(i,j+½) = ¼ dx [ a_z(i−1,j) + a_z(i−1,j+1) + a_z(i,j) + a_z(i,j+1) ]
Lz(i+½,j) = ¼ dz [ a_x(i,j−1) + a_x(i,j) + a_x(i+1,j−1) + a_x(i+1,j) ]
```

或者**直接用方案 B**（窗口本来就对）。G-3 的 45° 平面测试必须把"修正前 vs 修正后"都算一遍，
记录半格误差大小——这能反过来验证本判断。

### 2.2【P1】PML 的机制不是"单独处理"，是"叠加修正"

- `FiniteDifferenceSolver::EvolveE` 分派后对 **整个 Efield MultiFab** 做 `EvolveECartesian`
  （`EvolveE.cpp:95-100`，MFIter 无 PML 排除），然后 `WarpXPushFieldsEM.cpp:1048` 的
  `EvolveEPML` 在 PML box 内**再修正**（基于 Yee 系数的分裂场修正）。
- 新内核若照抄 box 逻辑，会在 PML 区产生"共形更新 + Yee 式 PML 修正"双重作用，PML 吸收失效甚至不稳定。
- 实现选项：(a) 内核外包一层"非 PML box 列表"（由 `m_do_pml_Lo/Hi` 构造的 box 集合取补集）；
  (b) `iMultiFab` 掩膜在核心里跳过 PML 单元。选 (a) 与现有 `EvolveEPML` 的 box 逻辑对齐。
- DESIGN2 的"EB 与 PML 相交则 abort"断言保留；再加一条：**新内核在 PML box 内不得生效**。

### 2.3【P2】方案 C 的 φ 符号

WarpX 的 `eb_implicit_function` **正值为固体**（`HANDOFF.md` 明确警告"写反会把外部全变导体"）。
方案 C 写 `Σ step·[φ>0]` 为开放长度——方向反了。开放 = `φ ≤ 0`（S1a' 用一个半开 cell 实测钉死）。
规则域上要保证采样积分给**精确** `dx/dz`，否则破坏 C-1 逐位。

### 2.4【P2】掩膜与 `A_y<=0` 不对齐

- 现状 `Ey` 掩膜：`ly==0` 即冻结；`ly` = "节点周围四条边任一 fully covered → 0"
  （`EmbeddedBoundaryInit.cpp:564-590`）。一个"三开放一覆盖"的节点会 `ly=0`、但 `A_y>0`。
- 新语义下二者冲突。建议：`Ey` 的更新判据与掩膜统一为 `A_y<=0`，并把 `MarkUpdateECellsECT`
  的 XZ `Ey` 分支改成同一判据；`Ex`/`Ez` 维持 `lx/lz==0`（与 `L<=0` 基本等价）。

### 2.5【P2】`F` 场（散度清洁 grad 项）不能丢

`EvolveECartesian` 在 Yee 更新后还有 `E += c²dt·grad(F)`（`EvolveE.cpp:203-228`）。若 ECT
打开了 div 清洁，新内核必须保留同一段。建议：新内核先照抄，并用"开关开 + F 启用"算例跑一次回归。

### 2.6【P2】F-2 自检空洞

`J=0` 时 `∂ρ/∂t+∇·J` 平凡为 0；`Σ(∇·B)` 检查的是 B 侧。**替代**：离散 Poynting 余量——
每步算 `ΔU = Σ(A_dual·E·ΔE 当量) + Σ(A_f·B·ΔB 当量)` 对照边界通量；无源真空腔里余量应 ≤ 机器零×步数
（几何精确时）。这才是"C 侧分子/分母不匹配"的直接探针。

### 2.7【P2】E-1 与稳定化/一阶几何冲突

`ΔU/U ≤ 1e-10` 只在"两侧几何同源且无稳定化"时现实（§3.2）。分档：

| 档 | 条件 | 目标 |
|---|---|---|
| E-1a | 方案 C（level set）+ 无稳定化 | ≤ 1e-10 |
| E-1b | 方案 A/B（一阶混合几何）+ 无稳定化 | 漂移随 dx 收敛（记录阶数） |
| E-1c | 打开 BCK | 有界、不发散即可；另测"抬升只发生在 f<f_c 的 cell" |

### 2.8【P3】小项

- §2.2 伪代码：`MaxwellSolverAlgo::ECT` → `ElectromagneticSolverAlgo::ECT`；`WARPX_DIM_XZ` 需 `#ifdef`。
- §5 3D 表：`Ex` 分子权重应为 `By×Ly、Bz×Lz`（对偶面 yz 的边是 y/z 方向）。
- D-10：2D 构建的 `IntVect` 只有 2 分量。
- §7-2"分子真实/分母稳定（等价于 B 侧的 Venl）"：不严格等价；Venl 是跨面借用记账，E 侧若真需要
  借用是"把抬升后的分母通量重新分摊给邻居分子"，是另一套记账。

---

## 3. 底层物理：对偶配对与离散能量

### 3.1 正确的 Hodge/FIT 表述（修订第一版）

- 半离散：`B -= dt·(⋆₂)^{-1} C E`，`E += c² dt·(⋆₁)^{-1} C̃ B`，其中 `⋆₂ = A_f`（B 面开放面积），
  `⋆₁ = A_dual`（E 对偶面开放面积），`C` 用主边开放长度 `l`，`C̃` 用对偶边开放长度 `L`。
- 能量守恒的充分条件（FIT/Weiland）：**`l`、`L`、`A_f`、`A_dual` 来自同一套 cut-cell 剖分**，
  使 `C̃ = C^T` 在带权内积下成立。第一版写"`A_dual·A_face = V²`、`C^T = −C`"都不成立，已弃用。
- 2D 的 `A_dual` 三种形态：`Ex`→线段 `Lz`、`Ez`→线段 `Lx`、`Ey`→面积 `A_y`。DESIGN2 已拆对。

### 3.2 现状与目标的配对差距

| 侧 | 分子几何 | 分母几何 | 同源性 |
|---|---|---|---|
| B（现状） | 主边 `lx/lz/ly`（`edge_cent`） | 主面 `face_areas`（`area_frac`/`vol_frac`） | 同一 EB 工厂，**同源** ✅ |
| E（新） | 对偶边 `Lx/Lz`（方案 A/B/C） | 对偶面 `Lz/Lx/A_y` | 与 B 侧**不同源**（除非方案 C 也用于 B 侧） ⚠️ |

- 方案 A/B 是"用主面/体积分数插值对偶量"的一阶近似 → 离散 Poynting 余量 O(dx)，能量漂移 O(dx)·步数。
- 要做到 E-1a（1e-10），需要 **B 侧的 `face_areas`/`edge_lengths` 也改用 level-set 同源积分**——
  这是"几何同源化"的独立工作项，本阶段可不做，但**必须写进文档**，否则会误以为 C 侧本身不守恒。

### 3.3 稳定化与守恒的取舍

- B 侧 BCK 已经"分母放大、分子真实"（`Rho = 环量/S_enlarged`），E 侧同款做法**对称**；
  DESIGN2 §2.5 判断正确。
- 抬分母必然破坏配对：BCK 开启后没有严格离散能量守恒，只有"有界 + 不发散 + 精度可控"（E-1c）。
- CFL：cut cell 上有效条件量级 `dt ≲ A_dual/(c·L_edge)`（2D 线段情形）；BCK 抬分母正是放宽它。
  Benkler 的闭式判据没有现成的 E 侧版本，**数值标定路线诚实**；建议扫描
  `f = 1e-1…1e-6 × 角度 0/15/30/45° × 分辨率 64/128`，把临界 `f_c` 与 `A_stab` 做成表。

---

## 4. 文献 / 教科书 / 开源 / 社区调研（增量）

### 4.1 论文

| 文献 | 要点 | 对本项目的意义 |
|---|---|---|
| Railton & Schneider 1999, "An analytical and numerical analysis of several locally conformal FDTD schemes", IEEE T-MTT 47(1):56–66 | 对多种局部共形格式做解析+数值分析，指出精度与稳定性不能从"更贴几何"自动推出 | **必须读**：本项目"BCK vs 共形 vs 阶梯"的三方案对比可以照它的框架做 |
| Dey & Mittra 1997（DOI 10.1109/75.622536） | 双侧共形 + 小胞元阈值 | 主线（沿用） |
| Xiao & Liu 2004 / 2008（DOI 10.1109/LMWC.2004.837384 / 10.1109/TAP.2008.916876） | ECT 扩胞 + 稳定性判据 | B 侧依据（沿用） |
| Benkler et al. 2006（DOI 10.1109/TAP.2006.875909） | BCK 主面判据 | E 侧 `A_stab` 无现成闭式，只能数值标定（§3.3） |
| Zagorodnov et al. 2007（JCP 225:1493–1507） | USC / 扩胞 | "借用"升级路径的备选 |
| Nieter et al. 2009（DOI 10.1016/j.jcp.2009.07.025） | VORPAL 实现 DM；dt 可任意缩小 | 小胞元代价的出处（沿用） |
| Kuo & Kuo 2010（[PDF，由 gprMax 报告引用](https://www.ee.nsysu.edu.tw/uploadPDF/PDF/P4162002884318.pdf)） | 共形 FDTD 小不规则胞元（cut 面积可无限缩小） | 米氏爆炸的独立文献佐证 |
| Esirkepov 2001（CPC 135:144–153）；Villasenor & Buneman 1992（CPC 69:306–316） | 任意形状因子的严格电荷守恒沉积 | J≠0（M2a-2）的必读；C 侧接口应预留开放面积 |
| Farjadpour 2006 / Kottke 2008 | subpixel（材料张量，非 PEC 几何） | 材料侧（M2b）；PEC 不能用它（§4.3 Meep 证据） |
| "Late-time growth is a serious defect of the time domain electromagnetic numerical methods"（[IEEE Xplore 检索到的作者主页摘录](https://ieeexplore.ieee.org/author/37347358500)） | 晚时增长是时域方法的公认缺陷 | 支持"晚时不稳定监控"建议（§6 S6）；正文未核到，仅作线索 |
| Stable non-orthogonal FDTD（[QML 论文页](https://eecs.qmul.ac.uk/~yang/documents/Eletters_FDTD.pdf)） | 非正交网格上的稳定格式 | "非 Yee 几何要单独做稳定性"的旁证 |

### 4.2 教科书

| 书 | 用途 |
|---|---|
| Taflove & Hagness, 3rd ed. 2005 | 共形 FDTD 章 + 晚时增长/稳定性的标准论述 |
| Gedney, *Introduction to the FDTD Method*, Morgan & Claypool 2011 | 共形格式稳定性分析（本项目的 `A_stab` 标定可参照其能量/谱框架） |
| Sullivan, *Electromagnetic Simulation Using the FDTD Method*, 2nd ed. 2013 | 共形方法的工程实现与常见坑 |
| Bossavit 1998 / LeVeque 2002 | Hodge 形式 / 守恒型离散（沿用） |

### 4.3 开源项目对照（本机源码级核查）

| 项目 | 核查事实（含出处） | 对本项目的结论 |
|---|---|---|
| **WarpX 上游** | E 侧 = 掩膜 Yee（`EvolveE.cpp:92-94`）；BCK 是"精度低于 ECT 的兜底"（`WarpXFaceExtensions.cpp:158-159` 注释） | C 侧在上游无参照实现（沿用第一版判断） |
| **AMReX** | `areafrac` = lo 面 MultiCutFab（`AMReX_EBDataCollection.cpp:60-64`）；`bndrycent`/`edgecent`/`volfrac` 可用；**对偶边无现成量** | 方案 A/B/C 都要自建；`bndrycent` 可作"方案 B+"（§6 S3） |
| **Meep** | `doc/docs/Subpixel_Smoothing.md` 原文：subpixel 只适用于**无频散无损介质**（"It only applies to frequency-independent, lossless dielectrics"）；由 **level-set** 解析算填充分数与法向；"directly discretizing a discontinuity…leads to 'stairstepped' interfaces"；**"fields directly on the interface are still at best first-order accurate"**；issue #1064：分散大时平滑可能伤稳定性 | ① 方案 C 与 Meep 同构（level set 解析）；② **PEC 不在 subpixel 范围**，须几何共形（支持本项目路线）；③ **界面上场仍一阶** → L3/Mie 判据应聚焦离界面量（反射/远场），管理预期 |
| **gprMax** | `testing/validation/impedance_surface/stability_report.md` 原文："The small-irregular-cell issue…is established for conformal FDTD, where continuously varying cut areas can shrink without bound"（引 Kuo & Kuo 2010）；"A retained voxel contributes one quarter of its volume to each incident electric edge and one half to each magnetic face"（**电对偶面积 ¼/½/¾ 构造先例**）；能量法引 Clemens & Weiland FIT（PIER 00080103）；PEC 是 `se == inf` 的**阶梯材料**（`gprMax/materials.py:248-260`），无 PEC 几何共形 | ① "电对偶面积按 ¼ 象限组装"与 `A_y` 的 ¼ 同构——**开源先例**；② gprMax 的共形只用于 SIBC 阻抗面，PEC 仍阶梯——印证"PEC 共形没有现成轮子" |
| VORPAL | Nieter 2009（闭源，论文可查） | 完整 DM 的唯一公开参考（沿用） |
| OpenEMS | 未核实其 PEC 处理 | 待核实（列入 §5 边缘项，不引用） |

### 4.4 社区 / issue

| 项 | 内容 | 用途 |
|---|---|---|
| WarpX [#7206](https://github.com/BLAST-WarpX/warpx/issues/7206) | 标题："ECT: B is never updated on faces stabilized by the BCK correction"（正文因网络未抓到，仅标题可证） | **警告**：E 侧做 BCK 时，`bck_stabilized` 的 E 点必须仍参与更新（用抬升后的分母），不要重蹈上游 B 侧覆辙 |
| Meep [#1064](https://github.com/NanoComp/meep/issues/1064) | subpixel 与稳定性的权衡（Subpixel_Smoothing.md:155 引用） | 平滑不保证稳——方案 C 上线后仍需稳定性扫描 |
| Clemens & Weiland, FIT（[PIER](https://www.jpier.org/PIER/pier.php?paper=00080103)） | 对称材料矩阵 → 能量守恒的严格框架 | §3.1 的正式出处 |

---

## 5. 疑问 / 存疑 / 质疑清单

> 每条给：依据、置信度（高/中/低）、不处理的后果。**这些不是"错误"，是必须边做边验证的开放项。**

| # | 疑问/存疑/质疑 | 依据 | 置信度 | 不处理的后果 |
|---|---|---|---|---|
| Q1 | Dey–Mittra 原文的 2D E 更新：`¼` 四象限是否**只**属于节点量？`Ex`/`Ez` 是否就是线段分母？ | 原文给 3D 对偶面；[UCalgary 实现说明](https://ucalgary.scholaris.ca/bitstreams/d5b1cb8f-0c73-4faf-9841-0c6171fc0932/download)称"E contours 不形变、按正常 FDTD 更新"；本文 §3.1 的推导 | 高（节点 ¼）、中（`Ex/Ez` 线段形式需 S1a' 手算复核） | 若 `Ex/Ez` 形式有误，TM 全部错 |
| Q2 | **质疑**：方案 A 的"位置二阶"只对 z/x 插值成立；对斜切面，volfrac/areafrac 与中心线开放分数的差是 O(dx)，方案 A/B 整体仍一阶 | 几何直观 + Meep 文档"界面上场一阶"的同源结论 | 高 | 把一阶方案当二阶用，收敛阶测试误判 |
| Q3 | BCK 的 `A_stab` 能否套 Benkler 主面判据换变量（面积→长度）？ | Benkler 判据针对主面；2D 对偶面是线段，形式不同 | 中 | 若硬套，临界 `f_c` 标定全错 |
| Q4 | **怀疑**：E 侧 BCK-only 能否根治米氏 80²/160² 爆炸？ | 上游 issue #7206 显示 BCK 面有"不被更新"的历史坑；BCK 是欠松弛 | 中 | S6' 若仍爆：直接试 DM 阈值丢弃 → 再升级借用/USC（§6 S10 已排好退路） |
| Q5 | 共形 E + Yee PML 修正（PML 区保持 Yee）在 EB 靠近 PML 时是否稳定？ | 无文献；DESIGN2 用 abort 规避；L3 算例 EB 距 PML 足够远 | 低 | 生产算例若 EB 近 PML 会静默出错；保留 abort |
| Q6 | 方案 C 在规则域的"精确常数"是否逐位可靠？ | parser 求值在规则域应全同号；但不同进程/向量化/采样边界的浮点一致性未验证 | 中 | C-1 逐位被浮点噪声破坏，误判索引错误 |
| Q7 | `Lx/Lz`（Ey 分子）与 B 侧 `edge_lengths` 是否必须同源？ | §3.2：不同源 → Poynting 余量 O(dx) | 高（已答，记录在案） | 能量漂移被误判为"C 侧不守恒" |
| Q8 | DESIGN2 §2.3.3"cut 但开放的点不再冻结"——现状对 `Ex/Ez` 已如此（`lx/lz==0` 才冻结），真正要改的只有 `Ey` 的 `ly` 判据 | `EmbeddedBoundaryInit.cpp:262-292` 代码事实 | 高 | 表述扩大化 → 多改 `MarkUpdateECellsECT`，引入无谓风险 |
| Q9 | 3D 的 `¼` 四象限平均对斜面的精度与 2D 相同吗？ | 3D 斜切时 face 分数与对偶面分数差异更大 | 中 | 3D 直接照搬 2D 的一阶方案可能更差；3D 优先方案 C |
| Q10 | `Fields.H` 加 3 个字段对旧 checkpoint/寄存器的影响 | 字段枚举变化影响所有 ECT 构建的 checkpoint 布局与字段索引 | 高（机制）；影响小 | 旧 checkpoint 读新二进制可能错位；需版本兼容处理或文档声明 |
| Q11 | E-1a（1e-10）在有 PML 的腔里能否达到？ | PML 是耗散层，无源腔能量不守恒 | 高 | 建议 E-1 用 **PEC 腔 + 解析腔模**（§6 S11） |
| Q12 | C 侧上线后 TM 基线 `r=−0.999361` 变好还是变坏？ | 变好：去掉 Ex/Ez 阶梯误差；变坏：一阶 `Lx/Lz` 引入新误差。两机制竞争 | 中 | 需 S4'/S5' 复测并解释方向 |
| Q13 | **质疑**：方案 C"640² 不到 1 秒" | parser 逐点求值 640²×2×16 ≈ 1.3e7 次 + ghost/MPI 边界重复，CPU 单核 0.1–1 s 量级；是估计不是承诺 | 高 | 把它当性能承诺写进计划会失真；按"≤ 数秒/层"表述 |
| Q14 | `getAreaFrac()[1]` 在 AMReX 2D 里对应哪个物理方向？ | 2D 的 dir=1 是第二维（z）；fork 的 `ComputeFaceAreas` 已按 `idim=2→amrex 1` 使用并过 L1 | 高 | 无风险，但 S1a' 应实测确认 |
| Q15 | 用 `bndrycent`（边界质心）算对偶段开放长度是否优于方案 A/B？ | AMReX 提供 `m_bndrycent`（`AMReX_EBDataCollection.cpp:55-58`）；直线界面可解析求交；未实测 | 中 | 错过一个几乎零成本的二阶改进（§6 S3） |
| Q16 | E 侧 BCK 的 `A_stab`/`L_stab` 临界值是否与 B 侧 `S_stab` 独立？ | 两侧几何量纲不同（面积 vs 长度）；Benkler 无 E 侧版 | 中 | 把 B 侧阈值直接搬来可能过强/过弱；分开标定 |

---

## 6. 建议清单（含理由、替代方案、代价）

| # | 建议 | 理由/依据 | 替代方案 | 代价 |
|---|---|---|---|---|
| S1 | S1a' 扩展为**几何真值表**：全开/全闭/左半/右半/下半/上半/对角 45°/135° 等 15 种单 cell 情形，手算开放长度写死在测试里 | lo/hi、φ 符号、四象限归属一次钉死（DESIGN2 自己的"不做后面全白费"） | 只测一种 | 0.5 天 |
| S2 | 默认方案 B；修正后的方案 A（§2.1 四面平均）做对照；G-3 输出两者的 45° 收敛曲线 | 方案 B 窗口正确；对比能量化"一阶"到底多大 | 只留一种 | 0.25 天 |
| S3 | 实现"方案 B+"：用 `bndrycent` 做直线界面的对偶段解析求交（Q15） | 若有效，以极小成本获得准二阶几何 | 跳过，留到方案 C | 0.5 天 |
| S4 | 判据聚焦**离界面量**（反射系数、远场、探针在界面外），不把界面上的逐点场当验收 | Meep 文档明说界面上场至多一阶（§4.3） | 继续用界面点 | 判据失真 |
| S5 | 做 Railton–Schneider 式三方案对比（阶梯 Yee / BCK / 共形 E+B）同一算例同一分辨率 | T-MTT 1999 的框架；直接回答"共形是否真更准更稳" | 只测新方案 | 0.5–1 天 |
| S6 | 晚时不稳定监控：每 1000 步记录 `max|E|`/`max|B|`/能量，任何 >1% 增长报警 | 晚时增长是时域共形方法的公认缺陷（§4.1） | 无监控 | 静默发散 |
| S7 | 把 issue #7206 教训写进 E 侧 BCK 实现注释：`bck_stabilized` 的 E 点**必须仍参与更新**（用抬升后的分母） | 上游 B 侧有历史坑 | 无 | 重蹈覆辙 |
| S8 | C 侧接口预留 `A_dual` 指针；J≠0 时 abort；将来先做 Esirkepov 主网格沉积 + 对偶面积修正，再做 cut-cell 守恒沉积 | Esirkepov 2001 / Villasenor–Buneman 1992；M2a-2 范围 | 现在就做守恒沉积 | 接口重写 |
| S9 | 新增 CTest：单 cell 几何真值表 + L2a 逐位 + `rotated_cube` 回归，2D 构建进 CI | 防回归；本 fork 目前无 C 侧测试 | 手工跑 | 0.5 天 |
| S10 | 预排升级路径：BCK-only → DM 阈值丢弃（`f<1e-3` 阶梯化）→ 借用/USC（Zagorodnov 2007） | Q4 的退路；每级都有文献依据 | 无退路 | 米氏卡死 |
| S11 | E-1 用 PEC 腔 + 解析腔模（无 PML） | Q11：PML 耗散破坏无源能量守恒 | 在 PML 腔里测 | 判据不可达 |

---

## 7. 全面性清单（第一版 + DESIGN2 都未覆盖，本次补齐）

| # | 遗漏项 | 建议 |
|---|---|---|
| 1 | 方案 A 半格偏移（§2.1） | 开工前改公式或用方案 B |
| 2 | PML box 排除机制（§2.2） | 内核外包非-PML box 循环；加断言 |
| 3 | 方案 C 的 φ 符号 | 用项目约定（正值=固体）；S1a' 实测 |
| 4 | `ly` 掩膜与 `A_y<=0` 对齐（§2.4） | `Ey` 判据统一为 `A_y<=0` |
| 5 | `F` 场 grad 项（§2.5） | 新内核保留该段 |
| 6 | F-2 自检空洞（§2.6） | 换成离散 Poynting 余量 |
| 7 | E-1 分档（§2.7） | E-1a/1b/1c 三档 |
| 8 | 几何同源性（§3.2） | 明确"E-1a 需 B 侧几何也同源"的工作项边界 |
| 9 | 2D `IntVect` 2 分量 | D-10 措辞修正 |
| 10 | 3D 表标签/半格 | 3D 开工前复刻 S1a'/G-3 |
| 11 | 重启后 `A_y` 除零保护 | 除 `Fields.H` 三项外，重启路径加"几何为 0 则重算或 abort" |
| 12 | TM 基线重记 | S4'/S5' 前先记旧基线 |
| 13 | 旧 checkpoint 兼容（Q10） | 字段枚举变更写进发布说明 |
| 14 | 方案 C 的确定性（Q6） | S1b' 在规则域验证 `max|L−dx|=0` 严格成立 |
| 15 | OpenEMS 的 PEC 处理方式 | 待核实；不引用 |

---

## 8. 开工条件与顺序建议

**开工前必改三处**：

1. 方案 A 修正（§2.1）或直接声明方案 B 为默认；
2. 内核的 PML box 排除设计（§2.2）；
3. F-2 换成 Poynting 余量、E-1 分档（§2.6/§2.7）。

**顺序**（DESIGN2 S0'–S11' 微调，含 §6 建议）：

```
S0'  提交基线 + 开分支（先做一次 git commit，fork 目前 0 提交）
S1a' 单 cell 真值表（§6 S1）：lo/hi 语义 + φ 符号 + "方案 A(原/修) vs 方案 B" 45° 解析对比
S1b' ComputeDualGeometry（方案 B 或修正 A），G-1/G-3
S2'  常数几何内核 → L2a 逐位（硬门禁）
S3'  真实几何 → L2a 仍逐位；Poynting 余量自检
S4'–S5'  L2b/L3 复测（往返时间；记录 TM 基线变化方向，Q12）
S6'  圆柱 C12 + 米氏 40/80/160（若仍爆 → §6 S10 退路）
S7'–S8'  稳定化（只抬分母）+ volfrac 扫描 + E-1 分档（PEC 腔模，§6 S11）
S9'  长跑 TM 1e4 步（漂移 ≤1e-4）
S10' 方案 C（level set，φ 符号按项目约定）→ 米氏 5 档收敛阶（阶数对照 §6 S4 的离界面量）
S11' 3D + 回归（先只开 2D 开关）
全程  晚时不稳定监控（§6 S6）
```

**硬门禁**：S2' 的 `max|d| = 0`。任何"看起来还行但差 1e-2"都是索引/符号错误，不许进入 S3'。

---

## 9. 参考文献

**本次新增（链接见 §4）**：Railton & Schneider 1999（T-MTT 47(1):56–66）；Gedney 2011；Sullivan 2013；
Kuo & Kuo 2010（gprMax 报告所引 PDF）；Esirkepov 2001（CPC 135:144–153）；Villasenor & Buneman 1992
（CPC 69:306–316）；Clemens & Weiland, PIER（FIT）；gprMax `stability_report.md`（本机
`<local-ref>/gprMax`）；Meep `Subpixel_Smoothing.md` 与 issue #1064（本机
`<local-ref>/meep`）；WarpX issue [#7206](https://github.com/BLAST-WarpX/warpx/issues/7206)（仅标题）；
[UCalgary Dey–Mittra 实现说明](https://ucalgary.scholaris.ca/bitstreams/d5b1cb8f-0c73-4faf-9841-0c6171fc0932/download)；
[ETH FIT 讲义](https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/154748/eth-46840-02.pdf)；
[arXiv:2601.01902](https://export.arxiv.org/pdf/2601.01902)；[HAL hal-03328215](https://hal.science/hal-03328215v1/document)；
[HAL tel-03820106](https://hal.science/tel-03820106v2/preview/2022LIMO0066.pdf)；
[AIP 1507(1):184](https://pubs.aip.org/aip/acp/article-pdf/1507/1/184/11481217/184_1_online.pdf)；
[IEEE Xplore 晚时增长线索](https://ieeexplore.ieee.org/author/37347358500)；
[QML stable non-orthogonal FDTD](https://eecs.qmul.ac.uk/~yang/documents/Eletters_FDTD.pdf)。

**沿用（DOI 已核）**：Dey & Mittra 1997（10.1109/75.622536）；Benkler et al. 2006
（10.1109/TAP.2006.875909）；Zagorodnov et al. 2007（JCP 225(2):1493–1507）；
Nieter et al. 2009（10.1016/j.jcp.2009.07.025）；Farjadpour et al. 2006（10.1364/OL.31.002972）；
Bossavit 1998；Taflove & Hagness 2005；LeVeque 2002；Birdsall & Langdon 1991。

---

## 附录 A　方案 A 修正速查

```
Lx(i,j+½) = ¼ dx [ a_z(i−1,j) + a_z(i−1,j+1) + a_z(i,j) + a_z(i,j+1) ]
Lz(i+½,j) = ¼ dz [ a_x(i,j−1) + a_x(i,j) + a_x(i+1,j−1) + a_x(i+1,j) ]
（或直接用方案 B：½ dx [v(i−1,j)+v(i,j)]、½ dz [v(i,j−1)+v(i,j)]）
```

## 附录 B　本次核对的关键行号

| 主题 | 位置 |
|---|---|
| 2D 交错 / ECT 分配 | `Source/WarpX.cpp:2459-2474`、`:2683-2750` |
| 2D E 模板 | `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp:167-199`（F 项 `:203-228`） |
| ECT 分派 | 同文件 `:92-94` |
| E 掩膜 XZ | `Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp:233-292`（`ly` 判据 `:564-590`） |
| `face_areas` XZ | 同文件 `:592-660`（`getAreaFrac` `:609/614`） |
| B 三分支 | `Source/FieldSolver/FiniteDifferenceSolver/EvolveB.cpp:285-395` |
| Rho 末尾重算 | `Source/FieldSolver/WarpXPushFieldsEM.cpp:1069-1090` |
| PML 修正调用 | 同文件 `:1048`（`EvolveEPML`） |
| checkpoint 字段 | `Source/Fields.H:133-135`、`:155-156` |
| AMReX areafrac lo 面 | `amrex-development/Src/EB/AMReX_EBDataCollection.cpp:60-64`（`bndrycent` `:55-58`） |
| gprMax SIBC 稳定性报告 | `<local-ref>/gprMax/testing/validation/impedance_surface/stability_report.md` |
| Meep subpixel 文档 | `<local-ref>/meep/doc/docs/Subpixel_Smoothing.md` |

*本报告与 `C_SIDE_DESIGN2.md` 配对使用；第一版的文献主线见 `C_SIDE_RESEARCH.md`；
B 侧审查见 `Docs/M2A_CODE_AUDIT.md`。*
</content>
