# B 调研报告：2D XZ 支持 TE 模（Bx/Bz）的共形 cut-cell 更新

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 范围：只做**源码阅读 + 文献/开源调研 + 设计问题**，**不实现**。
> 代码：`~/WarpX_fluid`（WarpX 26.10 分支）；参照 `~/WarpX_changed`（只读，逐字节一致）。
> 日期：2026-10-31

---

## 0. 结论摘要（先看这五条）

1. **2D XZ 的 ECT 是"面环"格式，只对 By 有意义**。源码事实：`face_areas[0]=face_areas[2]=0`、
   `edge_lengths[1]=0`（`EmbeddedBoundaryInit.cpp:533-537, 464-467`）；只有 `face_areas[1]`（y 面 = 元胞）
   非零。因此 `EvolveBCartesianECT` 里 Bx/Bz 因 `S<=0` 被直接跳过。
2. **2D 的交错是"E/B 共位"的对：Ex↔Bz 在 (i+½,j)，Ez↔Bx 在 (i,j+½)，Ey 在节点 (i,j)，By 在元胞中心 (i+½,j+½)**
   （`WarpX.cpp:2467-2477`）。所以 Bx/Bz 不是"面量"，而是"边量"。
3. **TE 的 Faraday 是 Ey 的 1D 导数，不是 Ez 的环量**：标准 Yee 在 2D 下是
   `Bx += dt·∂Ey/∂z`、`Bz += −dt·∂Ey/∂x`（`EvolveB.cpp:174,188`）。3D 的 `Rhox=(Ey·ly−…+Ez·lz−…)/Sx`
   在 2D 退化为 0/0（ly=0、Sx=0），**不能**直接复用。
4. **所以"补 2D 的 Rhox/Rhoz"不是照抄 3D 公式**，而是要设计一套**边/节点型共形离散**：
   用"cut 边长度"做 Bx/Bz 的 1D 共形 Faraday，用"节点对偶面积"做 Ey 的共形 Ampere。
5. 该工作**可以不改 WarpX 既有源码**（在 ECT 分支外新增一个 2D-TE 共形路径），但会与现有
   `EvolveB/EvolveECTRho` 的 TM 路径并存，需要设计好切换。**是否值得做，取决于目标**（见 §7 问题）。

---

## 1. 源码事实

### 1.1 2D XZ 场交错（`Source/WarpX.cpp:2467-2477`）

`IntVect` 分量顺序为 (x, z)；`1`=NODAL（节点/整数），`0`=CELL（元胞中心/半整数）：

| 分量 | nodal flag | 位置 | 角色 |
|---|---|---|---|
| Ex | (0,1) | (i+½, j) | x 边 |
| Ez | (1,0) | (i, j+½) | z 边 |
| Ey | (1,1) | (i, j) | 节点（出平面） |
| Bx | (1,0) | (i, j+½) | z 边（与 Ez 共位） |
| Bz | (0,1) | (i+½, j) | x 边（与 Ex 共位） |
| By | (0,0) | (i+½, j+½) | 元胞中心（出平面） |

要点：**Bx 与 Ez 同一位置，Bz 与 Ex 同一位置**；出平面的 Ey 在节点、By 在中心。

### 1.2 几何量在 2D 的处理

- `ComputeFaceAreas`（`EmbeddedBoundaryInit.cpp:513-566`）：
  - `idim==0 || idim==2` → `face_areas[idim] = 0`（2D 中 x/z 面是线段，零面积）；
  - `idim==1` → `face_areas[1] = getVolFrac()`（2D 体积分数即面积分数），再 `ScaleAreas` 乘 `dx·dz`。
- `ComputeEdgeLengths`（`:450-510`）：
  - `idim==1` → `edge_lengths[1] = 0`（2D 无 y 边）；
  - `idim==0/2` → 用 `getEdgeCent()` 得切边长度分数，再 `ScaleEdges` 乘 `dx/dz`。

### 1.3 为什么 Rhox/Rhoz 在 XZ 被跳过

`EvolveRhoCartesianECT`（`EvolveECTRho.cpp:111-146`）：`Rhox`/`Rhoz` 的赋值被 `#ifndef WARPX_DIM_XZ`
包住；且 `EvolveBCartesianECT`（`EvolveB.cpp:287`）对 `S<=0` 直接 `return`。两者一致地说明：
**2D 下 Bx/Bz 没有可用的面面积，ECT 面环格式不适用。**

---

## 2. 2D TE 的正确离散

### 2.1 方程组（Ey 出平面、Bx/Bz 入平面）

```
∂Ey/∂t = c²( ∂Bx/∂z − ∂Bz/∂x ) − Jy/ε0        (Ampere, Ey 在节点)
∂Bx/∂t = ∂Ey/∂z                                (Faraday, Bx 在 z 边)
∂Bz/∂t = −∂Ey/∂x                               (Faraday, Bz 在 x 边)
```

### 2.2 标准 Yee（规则网格，`EvolveB.cpp:174,188`）

```
Bx(i,j+½) += dt · [Ey(i,j+1) − Ey(i,j)] / dz
Bz(i+½,j) += −dt · [Ey(i+1,j) − Ey(i,j)] / dx
```

### 2.3 共形（cut cell）版本需要的新几何量

| 量 | 用途 | 现有量是否够 |
|---|---|---|
| z 边的 cut 开放长度（在 Bx 所在 (i,j+½)） | Bx 的 ∂Ey/∂z 分母/分子 | `edge_lengths[2]`（z 边）可提供，但取值位置需对齐 |
| x 边的 cut 开放长度 | Bz 的 ∂Ey/∂x | `edge_lengths[0]` 可提供 |
| 节点 (i,j) 的**对偶面积**（4 个相邻元胞开放面积之和） | Ey 的 Ampere 分母 | **没有现成量**，需由 `vol_frac` 构造 |
| PEC 上 Ey 的边界处理 | 切向 E=0 | 由 level set / 面质心判定，需新增 |

关键：**TM 用的是"面面积 + 边环量"，TE 用的是"边长度 + 节点对偶面积"**——两套几何量不同。

---

## 3. 为什么"补 Rhox/Rhoz"不成立

3D 的 `Rhox = (Ey·ly − Ey·ly + Ez·lz − Ez·lz)/Sx` 依赖：
- `ly` = y 边长度（2D 中恒为 0）；
- `Sx` = x 面面积（2D 中恒为 0）。

即便强行令 `Sx = dz`（把 x 面当作线段），公式里的 `Ez·lz` 项也与正确的 2D TE Faraday（用 **Ey** 的
z 导数）不符——3D 公式里 Bx 的环量由 Ez 构成，而 2D TE 里 Bx 的变化由 Ey 构成。**物理量张冠李戴**。
因此"补 Rhox/Rhoz"不能按 3D 公式平移，必须重写为边型 Faraday。

---

## 4. 三条可选实现路线

| 路线 | 做法 | 成本 | 风险/代价 |
|---|---|---|---|
| **B1 边/节点型共形** | 新增 2D-TE 共形 Faraday（cut 边长）+ 节点共形 Ampere（对偶面积），与 TM 的 ECT 并列 | 中高 | 需新几何量（节点对偶面积、PEC 上 Ey 掩膜）；与 ECT 并存需切换逻辑 |
| **B2 3D 虚拟层** | 把 2D 问题当成 1 层 3D（y 方向 1 格、周期）跑 3D ECT | 低 | 计算量 ~1 倍 3D；与 2D 构建不一致；仅作验证权宜 |
| **B3 直接用 3D 构建** | 2D 验证需要 TE 时改跑 3D（ECT 3D 无此限制） | 低 | 非 2D 能力；仅绕过 |

> 说明：B1 才是"真正让 2D 支持 TE"的路径；B2/B3 是绕过。

---

## 5. 稳定性（小胞元）

- TM 用小 cut **面面积** → dt 受面积分数限制，ECT 用面积扩展稳定。
- TE 用小 cut **边长度**（Bx/Bz 的 Faraday）与小 **节点对偶面积**（Ey 的 Ampere）→ 两类小量。
  类比 ECT，需要对"边长度"和"节点对偶面积"分别做扩展/借用。
- 这会引入第二套 `S_stab`/borrowing（边版与节点版），工作量不小，且需与 TM 的
  `ComputeSStab`/`ComputeFaceExtensions` 解耦或统一。

---

## 6. 验证算例（若实现）

1. **TE 法向入射平 PEC**：`Ey` 入射，测 `r = −1`（与 A 的 TM 基线同构）。
2. **TE 斜入射**：测反射角（2α）与能量。
3. **与 3D TE 对照**：同一物理，2D(XZ) vs 3D 薄层，逐点对比。
4. **收敛**：r 的误差随分辨率下降（对照 A 的 0.071%→0.032%→0.021%）。

---

## 7. 待决策问题（12 个，不同角度）

**目标与优先级**
1. B 的目标是"2D 也能跑 TE 以方便验证"，还是"2D 生产级 TE 能力"？（前者可用 B2/B3 绕过，后者才需 B1）
2. B 与 C（E 场共形 Ampere）的优先级：先 B 还是先 C？两者都涉及 E 场几何，是否应合并设计？
3. 是否接受"2D 只保证 TM、TE 一律走 3D"作为长期策略（即不做 B1）？

**实现与代码边界**
4. 是否允许修改 WarpX 既有源码（`EvolveB/EvolveECTRho/EmbeddedBoundaryInit`）？还是只在新模块
   `Source/MomentsSolver/` 内新增 2D-TE 路径（不改原码）？
5. 若新增 2D-TE 路径，它与现有 ECT(TM) 如何在运行时切换（按偏振？按 `algo` 新选项？）？
6. 是否接受"节点型 Ey"的共形 Ampere（与现有 cell-centered/edge-centered 更新风格不同）？

**几何与稳定化**
7. 节点对偶面积如何构造（4 元胞 `vol_frac` 求和？还是 `getVolFrac` 移位？）——是否已有 AMReX 原语？
8. TE 的小 cut 边长与小对偶面积，稳定化用"扩展/借用"还是"局部亚循环"？与 TM 的 ECT 是否统一？

**维度与范围**
9. 只做 2D XZ，还是也要 RZ（柱坐标）？RZ 的 TE 语义不同（m 模、轴奇点）。
10. 是否要覆盖 AMR（多层 TE）？还是单层先行？

**验收与成本**
11. TE 的验收标准（r=−1 精度目标、收敛阶）？是否要求与 3D 逐点一致？
12. B1 的工作量与"直接用 3D"的代价如何权衡——是否值得为"2D TE"单独投入？

---

## 8. 文献与证据索引

| 主题 | 来源 |
|---|---|
| 共形 FDTD 通用形式（面积比 A_f/A、边长比 ℓ_f/ℓ、PEC 切向 E=0） | Dey & Mittra 1997（DOI 10.1109/75.622536）；FDTD Basics（USPAS 讲义） |
| ECT（enlarged cell） | Xiao & Liu 2004（DOI 10.1109/LMWC.2004.837384）、2008（DOI 10.1109/TAP.2008.916876） |
| 2D ECT（弹性波类比，方法可迁移） | Yu & Mittra 2000；2D/3D ECT for elastic（GJI 2015/2016） |
| 2D 场交错 | `WarpX.cpp:2467-2477` |
| 2D 几何量 | `EmbeddedBoundaryInit.cpp:450-622` |
| ECT 面环与跳过 | `EvolveB.cpp:172-192, 287`；`EvolveECTRho.cpp:111-146` |
| 标准 Yee TE | `EvolveB.cpp:174,188` |
| WarpX 2D ECT 既有验证 | `Examples/Tests/embedded_boundary_rotated_cube`（仅 TM，checksum 已对齐上游） |

---

## 9. 结论

- **"补 Rhox/Rhoz" 的本质不是照抄 3D 面环，而是为 2D TE 设计一套边/节点型共形格式**；
  3D 面环公式在 2D 退化失效，物理量也不匹配。
- 是否做 B1，取决于"2D 是否真的需要生产级 TE"。若只是验证方便，4 的 B3（用 3D 跑 TE）成本最低。
- 建议先回答 §7 的 1–3 与 4–5（目标与代码边界），再决定是否投入 B1。

---

## 附录（2026-10-31 决策后修正）：per-unit-y 约化使 ECT 框架可直接扩展到 2D TE

> 用户决策（12 问答复）：2D/3D 都要生产级 TE；B+C 合并设计；可改 fork；稳定化用扩展/借用；
> 不做 RZ；覆盖 AMR；2D TE 须与 3D 逐点一致。据此修正 §3/§4 的悲观结论。

**修正**：§3 说"3D 面环公式在 2D 退化为 0/0、物理量张冠李戴"——那是对**当前代码几何取值**
（`face_areas[0]=face_areas[2]=0`、`edge_lengths[1]=0`）而言。若把 2D 几何按 **per-unit-y（单位 y 长度）**
解释，3D ECT 公式可以正确退化：

| 量 | 当前 2D 取值 | per-unit-y 正确取值 | 来源 |
|---|---|---|---|
| `face_areas[0]`（x 面，Bx 所在） | 0 | x-法向**线段**开放长度 = `getAreaFrac()[0] × dz` | AMReX EB2 2D 的 area_frac 即边长分数，位置恰为 Bx (i,j+½) |
| `face_areas[2]`（z 面，Bz 所在） | 0 | z-法向线段开放长度 = `getAreaFrac()[1] × dx` | 同上，位置恰为 Bz (i+½,j) |
| `edge_lengths[1]`（y 边，Ey 所在） | 0 | 节点开放度（流体=1，覆盖=0） | 节点 = 3D y 棱的 2D 退化 |

**验证约化**：3D `Rhox=(Ey·ly−Ey·ly+Ez·lz−Ez·lz)/Sx` 中，Ez 项因 y-不变性在两个 y 边界相消；
Ey 项（ly=1）除以 `Sx=dz`（或 cut 开放长度）正好给出 `∂Ey/∂z`——与标准 Yee 的
`Bx += dt·UpwardDz(Ey)`（`EvolveB.cpp:174`）一致。Rhoz 同理给出 `−∂Ey/∂x`。

**另需修复（独立于几何）**：`EvolveB.cpp:261` 的 ECT 循环 `for (idim=0; idim<AMREX_SPACEDIM; ++idim)`
在 2D 只到 idim=1（Bx、By），**Bz（idim=2）永不更新**——即使补了 Rhoz 也不会被消费；须扩为 `idim<3`。

**结论更新**：B1 改为 **B1'（在 ECT 框架内做 per-unit-y 约化）**，与 C（E 侧对偶几何共形 Ampere）
天然合并为同一套"对偶几何"设计（见 `cutcell-research.md` §8 修订）。§4 的 B2/B3 降级为
过渡手段。
