# C 侧设计与调研审查报告

> 审查对象：
> - `Docs/research_notes/C_SIDE_DESIGN.md`（198 行，2026-10-31 13:52）
> - `Docs/research_notes/C_SIDE_RESEARCH.md`（376 行，2026-10-31 15:10）
>
> 对照：WarpX fork `~/WarpX_fluid` 源码、上游 `~/WarpX_changed`、`HANDOFF.md`、`L3_REFLECTION_CRITERION.md`、
> `CASES_2026-10-31.md`；文献：Dey–Mittra、Benkler、Zagorodnov、Nieter、Farjadpour、Bossavit、Taflove、
> LeVeque、Birdsall–Langdon；开源：WarpX/AMReX、Meep、Gekyll（后者与 C 侧无关）。
>
> 审查方式：逐节读文档 → 回到源码核对 2D 交错与差分模板 → 重新推导 2D 对偶几何 → 与文献定义对照 →
> 检查基线数字与文档一致性。C 侧尚未实现，因此本报告审查的是**设计/调研的正确性与完备性**，不是运行结果。

---

## 0. 结论摘要

| 维度 | 结论 |
|---|---|
| 总体方向 | **正确**：E 侧仍是掩膜 Yee，单侧共形破坏 E/B 一致性；C 侧（对偶共形 Ampere）是必要修复。 |
| 文献框架 | **基本正确**：Dey–Mittra 双侧共形、BCK/USC/ECT、Meep subpixel、Hodge 复形，主线与引用大体成立。 |
| **核心几何** | **有系统性错误**：`Ex`/`Ez` 的 `A_dual=1`、切割记在"对偶边长"上，恰好把**分子/分母对调**；按文档公式实现，规则域退化会差一个 `dz`/`dx`，C-1 不可能通过。 |
| **Ey 内核** | **公式不对**：文档 §4.1 对相邻两个 `Bx` 取平均再 `/2`，不是 Yee 的"边中点值 × 开放长度"；规则域不能逐位退化。 |
| 存储设计 | **不完整**：`dual_face_area` 只有 (1,1) 一个节点数组，只够 `Ey`；`Ex`/`Ez` 的分母是各自的对偶边长，存储/索引未定义。 |
| 稳定化 | **未闭合**：`A_stab`/`L_stab` 是启发式，没有从 BCK 稳定性判据推导；E 侧不能照搬 B 侧的 `Venl`。 |
| 文献表述 | 若干数学表述需订正：`A_dual·A_face=V²`、`C^T=−C`、能量式缺 `ε0/μ0`、`dt∝√A_f` 的维度前提。 |
| 文档一致性 | `C_SIDE_DESIGN` 头部"L3 已通过"与 `C_SIDE_RESEARCH`/`HANDOFF` 的"30° 斜率 FAIL"冲突；L2b 指标混用绝对时间/往返时间。 |
| 完备性 | 缺：对偶几何的 AMReX 查询方案、掩膜新语义、PML/重启/AMR/2D-3D 的 C 侧版本、J 项接口、几何单测。 |

**一句话**：**调研（为什么做、文献依据、实施流程）质量高；设计（具体公式/索引/存储）在 Ex/Ez 对偶几何与 Ey 内核两处有关键错误，必须先订正再写代码，否则会在 C-1 上反复调试。**

---

## 1. 审查对象与方法

- **文档**：`C_SIDE_DESIGN.md`（目标/判据、几何量、实现清单、Ey 公式、稳定化、AMR、2D/3D、工作量）；
  `C_SIDE_RESEARCH.md`（问题清单、文献、单侧失稳分析、框架、S0–S10、风险、基线表）。
- **源码**：`EvolveE.cpp` 的 2D 差分模板与分派；`WarpX.cpp` 的 2D 交错标志与 ECT 分配；
  `InitializeEBGridData` 的 `maxLevel()` 限制；`EvolveB.cpp` 三分支；`EvolveECTRho.cpp`；`Fields.H`。
- **文献/开源**：核对 Dey–Mittra 的 E/H 双侧定义、BCK 稳定性判据、USC/ECT、VORPAL 实现、Meep subpixel、
  AMReX EB 数据（`volfrac`/`areafrac`/`edge_cent`）的语义。
- **未做**：没有实现 C 侧，也没有新增数值实验；基线数字取自 `HANDOFF.md`/`L3_REFLECTION_CRITERION.md`，
  仅做一致性核对。

---

## 2. 文档主张 ↔ 源码/推导核对

| # | 文档主张（位置） | 源码/推导事实 | 判定 |
|---|---|---|---|
| 1 | E 侧仍是掩膜 Yee，ECT 与 Yee 共用 `EvolveECartesian`（设计 §1；调研 P-C1） | `EvolveE.cpp:92` 的 `Yee || ECT` 分支；`:167-199` 三个分量模板 | ✅ |
| 2 | `m_eb_update_E` 只决定"是否更新"，不改变差分权重（设计 §1） | `:171/:182/:193` 的 `if (update_...==0) return;` | ✅ |
| 3 | 规则域上 ECT 与 Yee 必须逐位一致（C-1） | L2a 已证 B 侧；E 侧若模板正确也应如此 | ✅ 目标合理 |
| 4 | `Ex`/`Ez` 的 `A_dual = 1`，切割体现在对偶边长（设计 §3.1/§3.3；调研 §4.1） | 2D 模板：`Ex += c²dt(-DownwardDz(By))`（`:167-175`），分母是 `dz`；`Ez += c²dt(DownwardDx(By))`（`:191-199`），分母是 `dx` | ❌ **分子/分母对调**（详见 §3） |
| 5 | `Ey` 的 `A_dual = ¼ dx dz Σ volfrac`（设计 §3.2） | 节点对偶面由四个 ¼ 胞组成；正确 | ✅ |
| 6 | `Ey` 内核：`(Bx(i,j)Lz + Bx(i−1,j)Lz)/2 − …`（设计 §4.1；调研 §4.6） | Yee 的 `Ey` 模板是 `DownwardDx(Bz)` 与 `DownwardDz(Bx)`（`:180-188`），即**单个边中点值 × 开放长度**；没有相邻两点平均 | ❌ **不是 Yee**，C-1 不可能为 0 |
| 7 | 设计 §3.3 的退化检查 `A_dual=1, Lz=dz → /dz` | 代数上 `B·L/A = B·dz/1 = B·dz`，不是 `/dz`；自相矛盾 | ❌ |
| 8 | 新增数组：`dual_face_area` (1,1) + `dual_edge_length` 3×（设计 §3.2/§4.2） | `Ex`/`Ez` 的分母是对偶边长（不同位置），`Ey` 的分母是对偶面积；一个 (1,1) 数组不够 | ❌ 不完整 |
| 9 | 稳定化 `A_stab = 0.5·max{相邻对偶面(L·Δ)}`、`L_stab = 0.5·min{邻接 cell 完整边长}`（设计 §5；调研 §4.4） | B 侧 `ComputeSStab` 有严格定义（`WarpXFaceExtensions.cpp:96-133`）；E 侧没有对应推导 | ⚠️ 启发式，未闭合 |
| 10 | `A_dual·A_face = V²`、`C^T = −C`（调研 §2.1） | 对立方 cell，x 边的对偶面 `dy dz`、主面 `dy dz`，乘积 `(dy dz)² ≠ V²`；`C`（边→面）与 `C^T`（面→边）不是同一算子 | ❌ 数学表述错误/不精确 |
| 11 | Dey–Mittra 是双侧共形；WarpX ECT 只做 Faraday 一半（调研 §2.2） | 与 `EvolveE.cpp` 事实一致；DM 原文确实同时改 E/H | ✅ |
| 12 | E 侧"四个 cell 各 1/4"（调研 §2.2/§4.1） | 只对 2D 的**节点**（`Ey`）成立；`Ex`/`Ez` 是边，对偶面是线段，不是四象限 | ⚠️ 以偏概全 |
| 13 | 单侧共形 → 能量不守恒/局部增益（调研 §3.1） | 方向正确；具体式是量级估计（作者也注明） | ✅（启发式） |
| 14 | AMR 几何只在 `maxLevel()` 初始化（设计 §6） | `InitializeEBGridData` 的 `if (lev == maxLevel())`（`WarpXInitData.cpp:1608-1633`） | ✅ |
| 15 | 2D/3D 逐点一致只在规则域有先验保证（设计 §7） | 正确；cut cell 下 per-unit-y 与 3D 面环量无等价性证明 | ✅ |
| 16 | L3 已通过（设计头部） | `HANDOFF.md`：45°/60° 幅值 PASS，**30° 斜率 FAIL**（−10.5%，TM 横向性 0.446） | ⚠️ 表述误导 |
| 17 | L2b 等效壁面偏 `+0.37 ns`（约差一层 `dz`）（调研 §0） | `+0.371` 是 64² 的**往返**误差；等效壁面偏移 ≈ `c·0.371ns/2 ≈ 0.055 m ≈ 0.22 dz`，不是一层 | ⚠️ 数字对、解释过强 |
| 18 | S0 "新分支 `feat/ect-conformal-E`；`git status` 干净"（调研 §5） | fork 的 `master` **还没有任何提交**（42 个 untracked）；不先提交无法开分支 | ⚠️ 流程前置条件不成立 |

---

## 3. 核心几何：2D Ampere 的正确对偶量

### 3.1 源码事实（2D XZ，`EvolveE.cpp:167-199`）

WarpX 的 2D 交错（`WarpX.cpp:2459-2474`，`IntVect` 中 1 = nodal）：

| 量 | 位置 | 说明 |
|---|---|---|
| `Ex` | `(i+½, j)` | x 边（点 x 向，z 节点） |
| `Ey` | `(i, j)` | 节点（x、z 均节点） |
| `Ez` | `(i, j+½)` | z 边 |
| `Bx` | `(i, j+½)` | z 边（x 节点、z 格心） |
| `By` | `(i+½, j+½)` | 格心 |
| `Bz` | `(i+½, j)` | x 边 |

2D 模板（`DownwardDy ≡ 0`）：

```
Ex(i+½,j) += c²dt·( −∂z By )                 // :167-175，分母 dz
Ey(i,j)   += c²dt·( −∂x Bz + ∂z Bx )         // :180-188，分母 dx / dz
Ez(i,j+½) += c²dt·( +∂x By )                 // :191-199，分母 dx
```

### 3.2 正确的对偶量

由 3D Ampere 的环量形式（Dey–Mittra / Yee 的对偶面）在 ∂y=0 下退化：

| E 分量 | Ampere | 分母 `A_dual`（对偶面面积/单位 y） | 分子里的 B（× 开放长度） |
|---|---|---|---|
| `Ex` | `−c²∂z By` | **`Lz`：x=i+½ 处、z∈[j−½,j+½] 段的开放长度**（规则 `dz`） | `By(i+½,j±½)` × `1`（y 向不切） |
| `Ez` | `+c²∂x By` | **`Lx`：z=j+½ 处、x∈[i−½,i+½] 段的开放长度**（规则 `dx`） | `By(i±½,j+½)` × `1` |
| `Ey` | `c²(∂z Bx − ∂x Bz)` | **`A_y`：节点 (i,j) 周围对偶胞的开放面积**（`¼dx dz Σvolfrac`） | `Bx(i,j±½)×Lx`、`Bz(i±½,j)×Lz` |

要点：

1. **`Ex`/`Ez` 的切割量就是"对偶边长"**，但它出现在**分母**（作为对偶面面积），不是乘以 B 的权重；
   文档把 `A_dual` 写成 1、把 `L` 写成 `dz`，正好对调。
2. `Ey` 的 `Lx`/`Lz` 才是**乘以 B 的权重**；它的分母是 `A_y`。
3. 三个量有漂亮的对应关系：`Ex` 的分母 = `Ey` 的右/左对偶边 `Lz`；`Ez` 的分母 = `Ey` 的上/下对偶边 `Lx`。
   因此只需存 `Lx`、`Lz`、`A_y` 三套数组（`Lx`/`Lz` 在节点对偶胞的边上，`A_y` 在节点）。

### 3.3 文档错在哪

| 错点 | 文档写法 | 正确写法 | 后果 |
|---|---|---|---|
| ① `Ex`/`Ez` 的 `A_dual` | `A_dual = 1`（设计 §3.1/§3.3；调研 §4.1） | `A_dual = Lz`（`Ex`）/ `Lx`（`Ez`），规则 `dz`/`dx` | 退化后差 `dz`/`dx`；C-1 不可能通过 |
| ② 切割的归属 | "切割体现在对偶边长"（设计 §3.1） | 对 `Ex`/`Ez`：切割量是**分母**；对 `Ey`：切割量是分子权重 + 分母面积 | 公式整体错位 |
| ③ 退化检查 | "`A=1, Lz=dz → /dz`"（设计 §3.3） | `Lz=dz, A=Lz → /dz` | 文档内部矛盾 |
| ④ `Ey` 内核 | 相邻两点 `Bx(i,j), Bx(i−1,j)` 平均再 `/2`（设计 §4.1；调研 §4.6） | 单点边中点值 × 开放长度；`/2` 只来自 ¼ 胞的**面积/长度**，不是对 B 取平均 | 规则域不是 Yee；C-1 不可能为 0 |
| ⑤ 存储 | 一个 `dual_face_area` (1,1) + 3× `dual_edge_length`（设计 §3.2/§4.2） | 需要 `A_y`（节点）、`Lx`（z=j+½ 的 x 段）、`Lz`（x=i+½ 的 z 段）三套；`Ex`/`Ez` 的分母就是 `Lz`/`Lx` | `Ex`/`Ez` 无数组可读 |
| ⑥ 索引标签 | `Lx=areafrac[0]·dz`、`Lz=areafrac[1]·dx`（设计 §3.2） | `Lx` 是**沿 x** 的长度（z 法向面的开放长度），`Lz` 是**沿 z** 的长度（x 法向面） | 名称互换，易写错 |

**关于 ④ 的代数反例**（规则域，取 `Lx=dx, Lz=dz, A=dx dz`）：

文档公式给出

```
[ (Bx(i,j)+Bx(i-1,j))/2 − (Bx(i,j-1)+Bx(i-1,j-1))/2 ] / dz
−(同样的 Bz 项) / dx
```

即"相邻两个 z 差分的平均"；Yee 要的是 `(Bx(i,j)−Bx(i,j−1))/dz` 与 `(Bz(i,j)−Bz(i−1,j))/dx`。
两者相差一个 `x` 方向的平均（O(dx²) 的平滑），在规则域上**不为零**，因此 C-1 的 `max|d|=0` 不成立。
文档 §3.3 自己写的是单差分，与 §4.1 的公式不一致——这是设计文档必须首先修掉的自相矛盾。

### 3.4 对 C-1/C-2/C-3 的影响

- **C-1（规则域逐位）**：按现公式，`Ex`/`Ez` 差 `dz`/`dx` 量级，`Ey` 差 O(dx²) 平滑；不可能 PASS。
  这会直接卡住实施流程的 S2/S3。
- **C-2（L2b 壁面）**：E 侧变化会改变壁面等效位置；在公式修正前，测得的基线没有意义。
- **C-3（L3 幅值/斜率）**：同理；C 侧对 30° 斜率的改善无法在错误几何上评估。
- **C-0（开关关闭逐位）**：不受影响（走旧路径），但一旦打开就进入错误公式。

### 3.5 缺失的关键设计：对偶几何怎么从 AMReX 取

文档假设 `area_frac` 能给出 `Lx`/`Lz`（设计 §3.2；调研 §4.1）。**这不成立**：

- `Lz` 是 `x=i+½`（x 向**格心**）处的 z 段长度；`Lx` 是 `z=j+½`（z 向格心）处的 x 段长度。
- AMReX 的 `getAreaFrac()` 是**主面**（x=i 或 z=j）上的分数；`getVolFrac()` 是**格心**体积分数；
  `getEdgeCent()` 是**主边**（x=i+½/z=j 或 x=i/z=j+½）的质心。三者都不直接给"格心处的对偶段长度"。
- 需要的 `Lx`/`Lz` 是**对偶边**：例如 `Lz` 在 `x=i+½, z∈[j−½,j+½]`，它由两个**半主面**拼成，
  不是单个 `area_frac`。

可选方案（文档未写）：

1. **一阶**：`Lz ≈ 0.5·(volfrac(i,j)+volfrac(i+1,j))·dz`（或取相邻两格的较小/较大值）；
2. **精确**：用 level set / cut-cell 平面与直线 `x=i+½` 求交，积分开放长度（二阶，推荐用于最终精度）；
3. **对照**：在 45° 斜面上与解析开放长度比，量化一阶/二阶误差。

这是 C 侧真正的实现难点之一；文档把它简化成"读 `area_frac`"，会在 S1 就卡住或引入 O(dx) 误差。

---

## 4. 稳定化、能量与文献

### 4.1 稳定化

- 文档（设计 §5；调研 §4.4）提出 `A_stab = 0.5·max{相邻对偶面的 L·Δ}`、`L_stab = 0.5·min{邻接 cell 完整边长}`。
  - `A_stab` 的形式与 B 侧 `ComputeSStab` 同构，方向合理；但**没有推导**，也没有说明阈值系数 0.5 的来源。
  - `L_stab` 是多余的：`Ex`/`Ez` 的切割量是 `A_dual`（即 `Lx`/`Lz`），不需要再对"对偶边长"做一次稳定化；
    `Ey` 的 `Lx`/`Lz` 是分子权重，若被抬升会改变通量守恒性，需谨慎。
- 文献上，Benkler 2006 的稳定性判据是针对**主面**的 BCK；E 侧的 `A_stab` 需要单独推导（或做
  von Neumann/能量法分析），不能直接照搬。
- 文档"先 BCK-only、再决定借用"的谨慎是对的；但应先给出 `A_stab` 的闭式或数值标定方法。

### 4.2 能量

- 调研 §3.1 的 `U = ½Σ V_c|E|² + ½Σ A_f|B|²` 缺 `ε0`/`μ0`，是量级式；作为论证可以，但不能当严格判据。
- "单侧共形破坏 `C` 与 `C^T` 的伴随性"方向正确；但 `C^T=−C` 的写法不严格：`C`（边→面）与 `C^T`（面→边）
  是不同维度的算子，离散伴随关系要通过 Hodge star 表达（`C^T = ⋆^{-1} C^T ⋆` 之类），不是简单取负。
- C-5（100 周期 ≤1e-10）是合理目标，但需要定义能量泛函（含 `ε0`/`μ0`）与边界通量。

### 4.3 文献/开源对照

| 主题 | 文档表述 | 核对 |
|---|---|---|
| Dey–Mittra 双侧共形 | E/H 都改；Ampere 用对偶开放面积 | ✅ 主线正确；但"四个 cell 各 1/4"只适用于 2D 节点（`Ey`），`Ex`/`Ez` 是线段 |
| DM 小胞元阈值 | `A_f/A < α` 丢弃；`dt ∝ √(A_f)` | ✅ 方向对；`dt ∝ √A_f` 是 2D 说法，3D 应按最小边长/体积表述 |
| BCK | 抬 `A_f` 到 `A_stab`；WarpX `S_stab` 来源 | ✅；但 E 侧 `A_stab` 未推导 |
| USC/ECT | 扩胞换稳定；Zagorodnov 2007 | ✅；ECT 名称来源（IEEE 4463918）未核，但不影响 |
| VORPAL/Nieter | 完整 DM 实现；dt 可任意小 | ✅；Nieter et al., JCP 228 (2009) 7902–7916，DOI 10.1016/j.jcp.2009.07.025 |
| Meep subpixel | 材料界面平均；PEC 不是其目标 | ✅ 大体正确；Meep 也能用等效介质处理金属，但"PEC 几何共形"确非 subpixel 路线 |
| Hodge/复形 | `⋆₁`/`⋆₂`、`C`/`C^T` | ⚠️ `A_dual·A_face=V²` 与 `C^T=−C` 需订正 |
| AMReX EB | `volfrac`/`area_frac`/`edge_cent` | ✅ 但双几何查询未解决（§3.5） |

---

## 5. 接缝与实现清单核对

文档列出的 7 处接缝（调研 §4.5）与源码基本对得上：

| # | 文档位置 | 源码事实 | 备注 |
|---|---|---|---|
| 1 | `Fields.H` 加字段 | `Fields.H:70-72`、`:133-135` 现有 `edge_lengths/face_areas/area_mod` | 需同步 checkpoint 列表（`Fields.H:155-156` 已含 `ECTRhofield/Venl`） |
| 2 | `WarpX.cpp` 分配 | `:2683-2750` 的 ECT 分配段 | 仅 2D + 开关；要加 `FillBoundary` |
| 3 | `ComputeDualGeometry` | 新函数 | 几何查询方案未定（§3.5） |
| 4 | `InitializeEBGridData` | `WarpXInitData.cpp:1608-1633`，仅 `maxLevel()` | C 侧单层可先用；AMR 另做 |
| 5 | `EvolveE.cpp` 分派 | `:92` 的 `Yee || ECT` | 需加 `&& conformal_E` 分支 |
| 6 | 新内核 | 无 | 公式需按 §3.2 重写 |
| 7 | 对偶 `A_stab` | `WarpXFaceExtensions.cpp:96-133` 旁 | 未推导 |

**文档未覆盖的接缝**：

- **掩膜新语义**：C 侧上线后 `eb_update_E` 应只负责"covered 点强制 E=0"，不再承担精度职责；
  文档说"放在 S6 之后调整"，但未给出具体判据（哪些点必须写 0、哪些点参与共形更新）。
- **PML**：E 侧共形只在物理域；PML 仍走 Yee。需要确认 EB 不进入 PML，并加 TE+PML 回归。
- **重启/checkpoint**：新数组要加入 `Fields.H` 的 checkpoint 列表并在 `FillBoundary` 后一致；文档未提。
- **J 项**：接口保留（设计 §4.6 注 4），但未定义对偶面上的 J 积分（面积/边长权重）。
- **AMR/2D-3D**：方案有，但都排在 C 侧之后；可接受。

---

## 6. 文档一致性问题

| 问题 | 位置 | 事实 | 建议 |
|---|---|---|---|
| "L3 已通过" vs "30° FAIL" | 设计头部 vs 调研 §0/HANDOFF | 幅值判据（沿束平均 `|r+1|=0.0124`）PASS；角度/斜率判据（−10.5%、TM 横向性 0.446）FAIL | 改成"L3 幅值 PASS、30° 斜率 FAIL" |
| 0.0124 vs 0.026–0.051 | 设计 C-3 vs HANDOFF | 0.0124 是 30° 幅值；0.026–0.051 是 45°/60° | 标明角度与判据 |
| L2b 指标混用 | 设计 C-2（21.700 ns, 0.018 ns）vs HANDOFF（往返 +0.371/+0.017/+0.076 ns） | 绝对时间 vs 往返时间；两者都可用，但不能混 | 统一用往返时间（对载波相位不敏感） |
| `+0.37 ns ≈ 一层 dz` | 调研 §0 | 64² 往返误差 0.371 ns ⇒ 壁面偏移 ≈0.22 dz；"一层 dz"是 Yee 的壁面差 | 改为"等效壁面偏 ~0.2 dz" |
| S0 开分支 | 调研 §5 | `master` 无提交（42 untracked） | 先做一次基线提交，再开分支 |
| 工作量 5–6 天 | 调研 §7 | 未计入几何查询（§3.5）与 C-1 调试 | 加 1–2 天缓冲 |

---

## 7. 问题分级

| 级别 | 问题 | 影响 | 修复 |
|---|---|---|---|
| **P0** | `Ex`/`Ez` 的 `A_dual=1`、切割归属反了（设计 §3.1/§3.3；调研 §4.1） | 退化差 `dz`/`dx`；C-1 不可能 PASS | 按 §3.2 改：分母 = `Lz`/`Lx`，B 乘 1 |
| **P0** | `Ey` 内核的相邻点平均 `/2`（设计 §4.1；调研 §4.6） | 规则域不是 Yee；C-1 不可能为 0 | 单点边中点值 × 开放长度；去掉平均 |
| **P0** | 设计 §3.3 的退化检查算术自相矛盾 | 会误导实现 | 用 `L=dz, A=L` 重写 |
| **P1** | 双几何查询未定义（`areafrac` 不是对偶边） | S1 卡住或引入 O(dx) 误差 | 定义 `volfrac` 一阶 / level-set 二阶方案，并用 45° 斜面标定 |
| **P1** | 存储只够 `Ey`（`dual_face_area` (1,1)） | `Ex`/`Ez` 无分母数组 | 存 `Lx`、`Lz`、`A_y` 三套（或按 E 分量各一套） |
| **P1** | `A_stab`/`L_stab` 无推导 | 小胞元稳定性无保证 | 从 BCK/能量法推导或做扫描标定 |
| **P2** | 数学表述 `A_dual·A_face=V²`、`C^T=−C`、能量缺常数、`dt∝√A_f` | 文档可信度/推导严谨性 | 订正表述 |
| **P2** | 掩膜/PML/重启/J/AMR 的 C 侧版本缺失 | 生产化缺口 | 补设计 |
| **P2** | L3/L2b 指标与表述不一致 | 判据不可复现 | 统一指标 |
| **P3** | 引用细节（IEEE 文档号、Meep PEC 表述） | 小 | 核对 |

---

## 8. 建议（修正后的设计）

### 8.1 先改公式（P0）

**规则域退化目标**（必须逐位）：

```
Ex(i+½,j) += −c²dt·[By(i+½,j+½) − By(i+½,j−½)] / Lz(i+½,j)
Ez(i,j+½) += +c²dt·[By(i+½,j+½) − By(i−½,j+½)] / Lx(i,j+½)
Ey(i,j)   +=  c²dt·[ Bx(i,j+½)·Lx(i,j+½) − Bx(i,j−½)·Lx(i,j−½)
                    − Bz(i+½,j)·Lz(i+½,j) + Bz(i−½,j)·Lz(i−½,j) ] / A_y(i,j)
```

- 规则域：`Lz=dz`、`Lx=dx`、`A_y=dx dz` → 分别退化为 `−ΔBy/dz`、`+ΔBy/dx`、`ΔBx/dz−ΔBz/dx`，与
  `EvolveE.cpp:167-199` 一致。
- `Ex`/`Ez` 的分母就是 `Lz`/`Lx`（对偶边长），**没有**额外的 B 权重；`Ey` 才同时用 `Lx`/`Lz` 和 `A_y`。
- 不要对相邻 B 值取平均；`Lx`/`Lz` 的 ½ 只来自对偶胞的象限划分（如果采用分段积分）。

### 8.2 定义双几何查询（P1）

- `A_y(i,j) = ¼ dx dz [volfrac(i−1,j−1)+volfrac(i,j−1)+volfrac(i−1,j)+volfrac(i,j)]`（已有，正确）。
- `Lx(i,j+½)`：`z=j+½` 处、`x∈[i−½,i+½]` 的开放长度；由两个半主面拼接。
  - 一阶：`0.5·(areafrac_z(i−1,j)+areafrac_z(i,j))·dx`（注意索引约定）；
  - 二阶：level set 与直线求交，积分开放长度。
- `Lz(i+½,j)`：`x=i+½` 处、`z∈[j−½,j+½]` 的开放长度；同理。
- **验证**：45° 斜面，与解析开放长度比；规则域 `max|Lx−dx|=0`、`max|Lz−dz|=0`、`max|A_y−dx dz|=0`。

### 8.3 稳定化（P1）

- 先做 `A_stab` 的推导/标定：在 `volfrac` 扫描（1e-1…1e-6）中找临界 `f`，反推系数；
- `Ex`/`Ez` 只对分母 `Lz`/`Lx` 做 BCK；`Ey` 对 `A_y` 做 BCK，`Lx`/`Lz` 不抬（否则破坏通量守恒）；
- 保留 `C-0` 开关；先 BCK-only，再评估是否需要借用。

### 8.4 实施顺序（修正版）

1. **S0'**：先提交当前 B 侧代码（基线），再开分支；记录 8 组判据。
2. **S1'**：`ComputeDualGeometry`（`A_y`、`Lx`、`Lz`），只打印 min/max；规则域严格常数，45° 斜面与解析比。
3. **S2'**：常数几何 + **修正后的**内核（§8.1），跑 L2a，要求 `max|d|=0`；若不 0，先查符号/索引。
4. **S3'**：接真实几何，重跑 L2a（仍 0）。
5. **S4'–S9'**：L2b/L3/圆柱/稳定性/长跑/米氏，沿用 `C_SIDE_RESEARCH` 的 S4–S9，但基线用往返时间。
6. **S10'**：3D + 回归；AMR 单独一步。

### 8.5 文档修订

- 改 `C_SIDE_DESIGN.md` §3.1/§3.3/§4.1/§4.2/§5；
- 改 `C_SIDE_RESEARCH.md` §2.1 的 `A_dual·A_face=V²`/`C^T=−C`、§3.3 的"分母=对偶边长(L·dz)"、
  §4.1 的 `Ex/Ez` 行、§4.4 的 `L_stab`；
- 统一 L3/L2b 指标表述；修正"L3 已通过"与"+0.37 ns≈一层 dz"。

---

## 9. 文献核对

| 文献 | 文档引用 | 核对 |
|---|---|---|
| Dey & Mittra 1997 | MGWL 7(9):273–275 | ✅ DOI 10.1109/75.622536；双侧共形、小胞元阈值 |
| Benkler et al. 2006 | TAP 54(6):1843–1849 | ✅ DOI 10.1109/TAP.2006.875909；BCK 稳定性判据 |
| Zagorodnov et al. 2007 | JCP 225:1493–1507 | ✅（USC/扩胞） |
| Nieter et al. 2009 | JCP 228:7902–7916 | ✅ DOI 10.1016/j.jcp.2009.07.025；VORPAL 的 DM 实现 |
| Farjadpour et al. 2006 | OL 31:2972–2974 | ✅ DOI 10.1364/OL.31.002972；Meep subpixel |
| Bossavit 1998 | 书 | ✅ 离散形式/Hodge；但文档的 `A_dual·A_face=V²` 不是书中的表述 |
| Taflove & Hagness 2005 | 书 | ✅ Yee/稳定性 |
| LeVeque 2002 | 书 | ✅ FV/守恒 |
| Birdsall & Langdon 1991 | 书 | ✅ 电荷守恒（PIC 电流沉积） |
| WarpX glossary | ECT/BCK | ⚠️ IEEE 文档号未核；不影响技术结论 |

---

## 附录 A　修正后的 2D 对偶量表（可直接照做）

| 量 | 位置 | 规则值 | 用途 |
|---|---|---|---|
| `A_y` | 节点 `(i,j)` | `dx·dz` | `Ey` 的分母 |
| `Lx` | `z=j+½` 的 x 段（x∈[i−½,i+½]） | `dx` | `Ey` 上/下边权重；`Ez` 的分母 |
| `Lz` | `x=i+½` 的 z 段（z∈[j−½,j+½]） | `dz` | `Ey` 左/右边权重；`Ex` 的分母 |

**退化检查**（规则域）：

```
Ex: Lz=dz, ΔBy/dz ✓
Ez: Lx=dx, ΔBy/dx ✓
Ey: Lx=dx, Lz=dz, A=dxdz → ΔBx/dz − ΔBz/dx ✓
```

## 附录 B　行号索引（当前工作树）

| 主题 | 位置 |
|---|---|
| 2D 交错标志 | `Source/WarpX.cpp:2459-2474` |
| 2D E 模板（Ex/Ey/Ez） | `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp:167-199` |
| ECT 分派 | 同文件 `:92` |
| E 掩膜 | 同文件 `:171/:182/:193` |
| 2D 场分配 | `Source/WarpX.cpp:2683-2750` |
| ECT 初始化（仅 maxLevel） | `Source/Initialization/WarpXInitData.cpp:1608-1633` |
| B 三分支 / `idim<3` | `Source/FieldSolver/FiniteDifferenceSolver/EvolveB.cpp:264`、`:285-395` |
| Rho（XZ） | `Source/FieldSolver/FiniteDifferenceSolver/EvolveECTRho.cpp:120`、`:147` |
| `ComputeSStab`（XZ） | `Source/EmbeddedBoundary/WarpXFaceExtensions.cpp:96-133` |
| checkpoint 字段 | `Source/Fields.H:155-156` |

## 附录 C　本次审查未做的事

- 没有实现 C 侧，也没有跑新的数值实验；C-1/C-2/C-3 的判据仍是文档/HANDOFF 的基线。
- 没有逐行推导 BCK 稳定性判据；`A_stab` 的闭式需要专门做。
- 没有验证 AMReX `getAreaFrac()` 在 2D 的逐索引语义；建议在 S1' 用一个单 cell/45° 斜面手算钉死。

---

*本报告只审查 C 侧设计与调研文档；B 侧实现的审查见 `Docs/M2A_CODE_AUDIT.md`。*
</content>
