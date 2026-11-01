# M2a 合并设计的源码级核查：发现清单

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 目的：把 `cutcell-research.md` §8（B+C 合并设计）与 WarpX 26.10 fork 的实际源码逐条对照，
> 记录"代码事实 / 与文档的差异 / 需要触碰的位置 / 待决问题"。
>
> **本文不做最终拍板**：只记录事实与可选项；决策点在 §7 单列。
>
> 核查范围（逐文件通读，行号为当前工作树）：
> `EvolveB.cpp`(583)、`EvolveECTRho.cpp`(162)、`EvolveE.cpp`(570)、
> `WarpXFaceExtensions.cpp`(1014)、`WarpXFaceInfoBox.H`(121)、
> `EmbeddedBoundaryInit.cpp`(624)、`EmbeddedBoundaryInit.H`(143)、
> `WarpX.cpp` 的 ECT 分配段、`WarpXInitData.cpp` 的 ECT 初始化段、
> `WarpXPushFieldsEM.cpp` 的调用顺序、`WarpXEvolve.cpp` 的推进顺序、
> `EvolveEPML.cpp`/`EvolveBPML.cpp` 的 ECT 分支、`CartesianYeeAlgorithm.H` 的 2D 差分约定。
>
> 实验证据：`/tmp/ect_test/` 的 7 组受控实验；用户基线
> `benchmarks/cases/c02_0_pec/diags/reducedfiles/`；根因报告 `ROOTCAUSE_ECT_2D.md`。
>
> 核查日期：2026-10-31

---

## 0. 摘要（12 条事实）

1. **ECT 的 2D 实现是 TM-only**：`Rhox/Rhoz` 在 XZ 下不编译；`face_areas[0]/[2]` 与
   `edge_lengths[1]` 被硬置零；B 更新循环只到 `AMREX_SPACEDIM`；扩展/借用/BCK 只对 `idim=1`。
2. **2D TE 的失败链**：Bx 冻结（Rho=0 且 S=0）、Bz 从不更新 → Ey 只对冻结的 ∂zBx 做时间积分
   → 线性增长、不传播；激光 TE 偏振因此"零场"。
3. **2D TE 需要触碰的位置是 10 处**（§4 表），不是文档 §8.2(a) 写的"扩循环 + 补 Rhox/Rhoz"两处。
4. **E 侧对偶面/边是"平移半格"的新数组**：B 侧的对偶面恰好是主网格面（可用 `getAreaFrac`），
   E 侧不是；2D 下 Ex/Ez 的切几何体现在**对偶边长**上，Ey 体现在**对偶面积**上（§5）。
5. **§8.2(a) 的 2D B 侧公式与交错一致**（Bx 是 z 边、Bz 是 x 边，对偶面是主 x/z 面）——
   这一点是对的；上一版评审中"方向换反"的怀疑不成立（§6）。
6. **`EvolveE` 对 ECT 用的是同一套掩膜 Yee 内核**（`m_fdtd_algo == Yee || ECT`），
   仓库里不存在 `EvolveECartesianECT`。
7. **`EvolveECTRho` 在 E 更新末尾重算**（`WarpXPushFieldsEM.cpp:1069-1090`），
   供下一次 B 更新使用；三分支 B 更新用 `Rho` 与 `Venl` 记账（`EvolveB.cpp:285-381`）。
8. **3D ECT 对 TE/TM 都正常**（受控实验：3D TE 传播正常；2D TM 正常，2D TE 异常）。
9. **2D 的 `ComputeSStab<0>`/`MarkExtensionCells` 的 dim=0 公式访问 `k+1`**（2D 越界），
   dim=2 公式访问 `j+1`（合法）；但目前 2D 根本不调用它们。
10. **BCK 在 2D 被 `#if !defined(WARPX_DIM_XZ)` 禁用**（`WarpXFaceExtensions.cpp:597-620`）。
11. **`shrink_borrowing` 只处理 `idim < AMREX_SPACEDIM`**（2D 只到 1），`borrowing[2]` 未收缩。
12. **几何只在 `maxLevel()` 初始化/扩展**（`WarpXInitData.cpp:1607-1633`），
    AMR 多层需要另做（与 §8.6 的"先单层"一致，但代码事实要写清）。

---

## 1. ECT 数据流：初始化与每步

### 1.1 初始化（`WarpXInitData.cpp:1607-1633`，只在 `lev == maxLevel()`）

```text
ComputeEdgeLengths(edge_lengths, eb_fact);  ScaleEdges(edge_lengths, CellSize)
ComputeFaceAreas (face_areas,  eb_fact);    ScaleAreas (face_areas,  CellSize)
MarkExtensionCells(cell_size, flag_info, flag_ext, B, face_areas, edge_lengths, area_mod)
ComputeFaceExtensions()          // one-way → eight-ways → BCK（并打印警告）
MarkUpdateECellsECT(m_eb_update_E, edge_lengths)
MarkUpdateBCellsECT(m_eb_update_B, face_areas, edge_lengths)
```

- 数组分配在 `WarpX.cpp:2683-2750`：`edge_lengths[idim]` 在 **E_idim 的节点/边交错**上，
  `face_areas/area_mod/Venl/ECTRhofield[idim]` 在 **B_idim 的交错**上；
- `face_areas` 的初值是 0，`area_mod` 在 `MarkExtensionCells` 里被置为 `S`
  （`EmbeddedBoundaryInit.cpp:403`）；
- 扩展链在 `WarpXFaceExtensions.cpp`：
  `ComputeFaceExtensions`（`:548`）→ `ComputeOneWayExtensions`（`:635`）→
  `ComputeEightWaysExtensions`（`:795`）→ `ApplyBCKCorrection`（`:168`，兜底）。

### 1.2 每个时间步（`WarpXEvolve.cpp:646-681`）

```text
EvolveB(0.5 dt)   → FillBoundaryB
EvolveE(dt)       → … → EvolveECTRho（在 EvolveE 末尾，E 更新并施加边界之后）
FillBoundaryE
EvolveB(0.5 dt)   → PML 阻尼
```

即：**Rho 由"刚更新完的 E"算出，供下一次 B 更新（以及同一步的后半）使用**；
B 的三分支更新在 `EvolveB.cpp:285-381`，`Venl` 每次调用先清零（`:246-248`）。

### 1.3 B 更新的三分支（`EvolveB.cpp`）

| 分支 | 条件 | 行为 |
|---|---|---|
| unstable | `flag == extended` | `Venl = Rho*S + Σ Rho_neigh*borrowed_area`；`rho_enl = Venl/S_mod`；`B -= dt*rho_enl`；并把 `rho_enl*borrowed_area` 原子加回 `Venl` 的邻居（`:285-362`）|
| stable（未借出） | `available` 或 `bck_stabilized` | `B -= dt*Rho`（`:370-374`）|
| stable（已借出） | `intruded` | `Venl += Rho*S_mod`；`B -= dt*Venl/S`（`:375-379`）|

- `S = face_areas`，`S_mod = area_mod`；`S <= 0` 直接返回（`:287`、`:366`）；
- 借用邻居索引由 `FaceInfoBox::uint8_to_inds` 给出二维偏移，再按 `idim` 映射到三维坐标
  （`:299-319`、`:331-351`）。

---

## 2. 2D XZ 的实际约定（读源码后确认）

### 2.1 尺寸与索引

- `WarpX::CellSize` 在 XZ 返回 `{dx, 1.0, dz}`（`WarpX.cpp:3272-3286`）：
  **3-矢量的分量 1 是 y（=1），分量 2 是 z（=dz）**；
- 数组的第二维 j 是 z（`CartesianYeeAlgorithm::UpwardDz` 在 XZ 用 `F(i,j+1,k)-F(i,j,k)`，
  `CartesianYeeAlgorithm.H` XZ 分支）；
- 2D 的 3 个分量交错（`WarpX.cpp:2469-2474`）：

| 量 | 2D 交错 | 位置 |
|---|---|---|
| Ex | (0,1) | x 边：`(i+1/2, j)` |
| Ey | (1,1) | 节点：`(i, j)` |
| Ez | (1,0) | z 边：`(i, j+1/2)` |
| Bx | (1,0) | **z 边**：`(i, j+1/2)` |
| By | (0,0) | 格心 |
| Bz | (0,1) | **x 边**：`(i+1/2, j)` |

- 这一交错在 `MarkUpdateBCellsECT` 的注释里被明确写出：
  "In XZ and RZ, Bx lives on a z-edge"（`EmbeddedBoundaryInit.cpp:329`）、
  "Bz lives on a x-edge"（`:342`），并用 `lz`/`lx` 作覆盖判据（`:330`、`:343`）。

### 2.2 2D 下被硬置零的量

| 量 | 现状 | 位置 |
|---|---|---|
| `edge_lengths[1]`（y 向边，Ey 用） | `setVal(0.)`，`ScaleEdges` 跳过 | `EmbeddedBoundaryInit.cpp:465`、`:581` |
| `face_areas[0]`、`face_areas[2]`（Bx、Bz 的对偶面） | `setVal(0.)`，`ScaleAreas` 跳过 | `:535`、`:608` |
| `flag_info_face[0]/[2]`、`flag_ext_face[0]/[2]` | 在 `MarkExtensionCells` 里置 0 并 `continue` | `:376-382` |
| `area_mod[0]/[2]` | 因上一条从未初始化 | `:395-403` |

### 2.3 2D 下只处理 By 的四处循环

| 位置 | 现状 |
|---|---|
| `WarpXFaceExtensions.cpp:42-44`（`CountExtFaces`） | `for idim=1; idim<AMREX_SPACEDIM` |
| `:654-656`（`ComputeOneWayExtensions`） | 同上 |
| `:814-816`（`ComputeEightWaysExtensions`） | 同上 |
| `:597-604`、`:613-620`（BCK） | idim 0/2 被 `#if !defined(WARPX_DIM_XZ)` 排除 |

### 2.4 3D 公式在 2D 的越界/错位点

| 位置 | 问题 |
|---|---|
| `ComputeSStab<0>`（`WarpXFaceExtensions.cpp:97-100`） | 用 `ly(i,j,k+1)`：2D 的 k 恒 0 → k+1 越界（当前 2D 不调用） |
| `MarkExtensionCells` idim=0（`EmbeddedBoundaryInit.cpp:405-407`） | 同上（`ly(i,j,k+1)`） |
| `GetNeigh/SetNeigh/GetNeighPtr` dim=0（`WarpXFaceExtensions.cpp:280-306`、`:331-360`、`:382-408`） | 用 `arr(i, j+i_n, k+j_n)`：k 是 3D-y，2D 下错位 |
| `EvolveB` 的 idim=0 借用映射（`EvolveB.cpp:299-302`、`:331-334`） | `kp = k + vec(1)`：同上 |
| `shrink_borrowing`（`WarpXFaceExtensions.cpp:246`） | 只到 `AMREX_SPACEDIM`，2D 下 `borrowing[2]` 不收缩 |

---

## 3. 2D TE 为什么失败（因果链 + 实验）

### 3.1 因果链

```text
face_areas[0] = face_areas[2] = 0
  → EvolveB 的 idim=0/2 在 `if (S <= 0) return;`（EvolveB.cpp:287/366）处直接返回
  → Bx 恒为初值；Bz 连循环都不进（idim < AMREX_SPACEDIM = 2）
Rhox = Rhoz = 0（EvolveECTRho.cpp:117-120、141-144 在 XZ 不编译）
  → 即使 S > 0，Bx 也只会做 `B -= dt*0`
EvolveE 仍按 ∂t Ey = c²(∂z Bx − ∂x Bz) 更新（EvolveE.cpp:184-187）
  → Bx 有 ∂z 梯度但冻结 → Ey 线性累积（不传播）
```

### 3.2 受控实验（`/tmp/ect_test/`，同一 Gaussian 脉冲，2D 128²/3D 64³）

| 实验 | 求解器 | 偏振 | 观察 |
|---|---|---|---|
| `run_te` | 2D ECT | TE（Ey+Bx） | Bx = 3.315e-3 **逐位不变**；Ey 9.879e5 → 1.553e7（400 步）；质心 −2.02 不动 |
| `run_tm` | 2D ECT | TM（Ex+By） | 质心 −2.02 → −0.27（t=5.897e-9 s，c·t=1.768 m）✓ |
| `run_yee` | 2D Yee | TE | 质心 −2.02 → −0.23（1.79 m）✓ |
| `r3d_te` | **3D ECT** | TE | 质心 −2.03 → +0.91（t=9.629e-9 s，c·t=2.887 m）✓ |
| `l_te_ect` | 2D ECT + 激光 | (0,1,0) | 观察区 Ex=Ey=Bx=By=**0** |
| `l_tm_ect` | 2D ECT + 激光 | (1,0,0) | Ex = 1.86e8 ✓ |
| `l_te_yee` | 2D Yee + 激光 | (0,1,0) | Ey = 1.86e8 ✓ |

### 3.3 用户基线数据（`c02_0_pec/diags/reducedfiles/`，401 点，0 → 5.90e-8 s）

| 量 | 值 |
|---|---|
| `fp_inc`（x=0, z=−5）Ey | 1.83e4 → 8.83e5，**r² = 1.0000**，斜率 1.466e13 V/m/s |
| `fp_inc` Bx | **−6.1385e-5 恒定**（等于 t=0 初值） |
| `fp_shadow`（z=+4） | 全部分量恒 0 |
| 解析预言（该算例自己的包络） | $c\,|\partial_z E_{y0}| = c\cdot\frac{4}{1.5}e^{-4}\cdot10^6 \approx 1.47\times10^{13}$ V/m/s |

实测与解析预言相差 < 0.3%，证明"B 冻结、E 只对冻结的 ∂zBx 做时间积分"。

---

## 4. 2D TE 需要触碰的位置（10 + 2 处）

> 这些是"若要 2D 全 TE"必须改的点；行号即当前工作树。

| # | 文件:行 | 现状 | 改动性质 |
|---|---|---|---|
| 1 | `EmbeddedBoundaryInit.cpp:376-382` | XZ 下 idim 0/2 置 0 并 continue | 去掉跳过 |
| 2 | `EmbeddedBoundaryInit.cpp:465` | `edge_lengths[1] = 0` | 改为节点开放度（每单位 y） |
| 3 | `EmbeddedBoundaryInit.cpp:535` | `face_areas[0]/[2] = 0` | 改为 `getAreaFrac(0)×dz`、`getAreaFrac(1)×dx` |
| 4 | `EmbeddedBoundaryInit.cpp:581`、`:608` | `ScaleEdges` 跳过 1；`ScaleAreas` 跳过 0/2 | 去掉跳过 |
| 5 | `EmbeddedBoundaryInit.cpp:405-419` | idim=0 用 `k+1`；idim=2 用 `j+1` | 2D 版 S_stab（dim=0 把 `k+1` 改成 `k`） |
| 6 | `WarpXFaceExtensions.cpp:44`、`:656`、`:816` | 循环 `idim=1..` | 扩到 0/2 |
| 7 | `WarpXFaceExtensions.cpp:97-121` | `ComputeSStab<0>/<2>` | 加 2D 分支（dim=0 的 k+1） |
| 8 | `WarpXFaceExtensions.cpp:280-306`、`:331-360`、`:382-408` | dim=0 用 3D-y 索引 | 加 XZ 分支 |
| 9 | `WarpXFaceExtensions.cpp:597-620` | idim 0/2 的 BCK 被 `#if` 禁用 | 2D 启用 |
| 10 | `EvolveB.cpp:261` | 循环到 `AMREX_SPACEDIM` | 扩到 3 |
| 11 | `EvolveB.cpp:299-302`、`:331-334` | idim=0 借用索引 `k+vec(1)` | 加 XZ 分支 |
| 12 | `EvolveECTRho.cpp:117-120`、`:141-144`、`:155` | Rhox/Rhoz 在 XZ 不编译 | 补 XZ 公式与 `ignore_unused` |

**另有 2 处与 Bz 新增相关、容易被忽略：**

- `WarpXFaceExtensions.cpp:246` 的 `shrink_borrowing` 只到 `AMREX_SPACEDIM`，2D 下 `borrowing[2]`
  不会被收缩（`init_borrowing` 在 `:219` 是 0..2）——若 Bz 走借用，必须补；
- `CountExtFaces` 的归约 `ReduceIntSum(&(sums(0)), AMREX_SPACEDIM)`（`:66`）只归约 2 项，
  若统计 0/2 需归约 3 项。

---

## 5. E 侧共形 Ampere 的几何需求（源码 + 2D/3D 交错推导）

### 5.1 3D（参考 Xiao–Liu 2008 的成对构造）

- E 在边上；Ampere 更新需要 **E 边的对偶面面积**（dual face）与 **对偶边长**（dual edges）；
- 对偶面 = 主网格的"对偶面"（位于主网格格心/边心），不是 `face_areas`（主 B 面）；
- 对偶边长 = 对偶面边界上各段的开放长度，也不是 `edge_lengths`（主 E 边）。

### 5.2 2D 逐分量（本报告的新发现）

| E 分量 | 位置 | Ampere 的切向微分 | 2D 下被 EB 切割的几何量 |
|---|---|---|---|
| Ex | x 边 `(i+1/2, j)` | $\partial_t E_x=-c^2\partial_z B_y$（$\partial_y=0$） | **对偶边长**：x=i+1/2 处 z 段的开放长度 |
| Ez | z 边 `(i, j+1/2)` | $\partial_t E_z= c^2\partial_x B_y$ | **对偶边长**：z=j+1/2 处 x 段的开放长度 |
| Ey | 节点 `(i,j)` | $\partial_t E_y=c^2(\partial_z B_x-\partial_x B_z)$ | **对偶面积**：以节点为中心的 x-z 对偶 cell 的开放面积（4 个相邻 cell 各贡献 1/4，**不是** 4 个面积直接相加） |

- 关键事实：2D 下 Ex/Ez 的对偶面（y-z / x-y 面）在 y 方向无限，**面积 = dy = 1，不被切割**；
  切割只体现在**对偶边长**上；Ey 的切割体现在**对偶面积**上。
- 因此 §8.2(b) 的公式 `E += dt·(c²∮B·dl_dual − J·A_dual/ε0)/A_dual` 在 2D 要按分量解释：
  Ex/Ez 的 $A_{\rm dual}=1$（每单位 y），$dl_{\rm dual}$ 是开放段长；Ey 的 $A_{\rm dual}$ 是 open area。
- 这些对偶量都在**半格偏移**的网格上（cell center / node），与主 `edge_lengths`/`face_areas`
  的采样点不同，必须新建数组（例如 `dual_edge_length[3]`、`dual_face_area`）。

### 5.3 J 项

- `EvolveE` 现在直接减 `μ0*j`（`EvolveE.cpp:176/187/198`），J 沉积在主 Yee 边/面上；
- 共形版需要 $\int_{\rm dual} J\cdot dA$：在 2D 下 Ex 的 Jx 与 Ex 同点，但双面对应的开放段
  可能被切；Ey 的 Jy 在节点，对偶 cell 的 open area 可能被切；
- `EvolveECTRho` 在 E 更新**之后**重算（`WarpXPushFieldsEM.cpp:1069-1090`），
  所以 B 侧看到的是共形 E 的环量；J 与 Rho 的一致性（∇·E−ρ）需要单独检验。

### 5.4 与现有掩膜的交互

- `MarkUpdateECellsECT` 用 `lx/lz` 给 Ex/Ez 打掩膜，用节点周围四条边给 Ey 打掩膜
  （`EmbeddedBoundaryInit.cpp:262-288`）；
- 换成共形 E 更新后，掩膜的语义（"被覆盖的 E 不更新"）仍可保留，但 Covered 区的 E 值
  需要定义（PEC 上 E_t=0），且与对偶面积为零的处理要一致。

---

## 6. 对 `cutcell-research.md` §8 的逐条核对

| 文档条目 | 与源码的对照 | 结论 |
|---|---|---|
| §8.0 "B 侧用 B 面面积 + E 边长；E 侧用对偶面面积 + 对偶边长" | 与源码一致；3D 如此 | ✅ 认识正确 |
| §8.2(a) `face_areas[0] = getAreaFrac()[0]×dz`、`face_areas[2] = getAreaFrac()[1]×dx` | 与 2D 交错一致（Bx 是 z 边、Bz 是 x 边；对偶面是主 x/z 面） | ✅ 公式正确；但**未写**"当前代码把这些置 0" |
| §8.2(a) `edge_lengths[1] = 节点开放度（1/0）` | 与"每单位 y"一致（`CellSize[1]=1`） | ✅ 方向正确；但**未写**当前是 `setVal(0)` |
| §8.2(a) 循环扩 `idim<3` | 必要但远不充分（§4 列了 12 处） | ⚠️ 不完整 |
| §8.2(b) "对偶面面积 ≠ face_areas；对偶 S_stab/borrowing/Venl 记账/BCK" | 方向正确；但未给存储设计、未区分"2D 下 Ex/Ez 的切几何在对偶边长" | ⚠️ 需细化 |
| §8.2(b) "2D TE：Ey 的节点对偶面积（4 邻胞开放面积和）" | 字面会大 4 倍；应为 Σ(1/4)·A_open（或精确积分） | ⚠️ 需订正 |
| §8.2(c) "两侧统一 ECT 式扩展，BCK 前核对 #7206" | 与代码一致 | ✅ |
| §8.3 流体 "EB vol_frac + SRD + 源项分裂；PEC 壁 ghost-cell + 镜像 Riemann" | 与 AMReX SRD/文献一致；源码侧未接 | ✅（实现待做） |
| §8.4.1 "2D TM 基线 r=−0.9993 作回归" | 换全共形 E 后该值会变 | ⚠️ 应记为"改动前基线" |
| §8.4.2 "2D TE 与 3D 逐点一致（硬判据）" | 需定义 3D 设置（1 层 y、周期、dy=1、同 dx/dz/dt）与容差 | ⚠️ 判据不可执行直到定义 |
| §8.4.4–8.4.7 | 与文献/实验一致 | ✅ |
| §8.5 子任务 M2a-1a/1b | 相对 §4 的 12 处，1a/1b 的范围偏小 | ⚠️ 需重划 |
| §8.6 "设计覆盖 AMR、先单层" | 与源码一致（几何只在 maxLevel 初始化/扩展） | ✅（多层待做） |

---

## 7. 待决问题（不拍板）

> 以下每条只列事实、可选项与代价，供你决定。

### D1　2D TE 的时间位置

- 事实：2D TE 要动 12 处，分布在 3 个文件；3D ECT 对 s/p 已经正常。
- 选项 A：3D 先做 E 侧（3D 的 B 侧已通过 s/p 验证），2D TE 作为后续独立里程碑；
- 选项 B：2D TE 先做，3D 后做；
- 选项 C：两者并行（需要两套验证基线）。
- 代价差异主要在：2D TE 的改动面 vs 3D 每次运行的机时；以及"2D TE 与 3D 一致"这条判据
  只有在 3D 侧可信时才成立。

### D2　E 侧对偶量的构造方式

- 选项 A：从 level set 精确积分对偶面/对偶边（二阶，需访问 EB 的 level set 或 cut cell 几何）；
- 选项 B：用体积/面积分数近似（如 Ey 的 4×1/4 求和、Ex 的对偶边长用相邻 cell 的分数平均）；
- 事实：WarpX 现有的 `getAreaFrac/getVolFrac` 只给主网格量，选项 A 需要新接口。

### D3　2D TE 的稳定化范围

- 选项 A：完整 ECT（借面积 + 八邻域 + BCK），即 §4 的全部 12 处；
- 选项 B：先只做"无借用 + BCK 兜底"（跳过 one-way/eight-ways），验证 TE 传播后再补借用；
- 事实：BCK 在 2D 当前被 `#if` 禁用（`WarpXFaceExtensions.cpp:597-620`），选项 B 也要改这一处。

### D4　J 项

- 选项 A：先做 J=0 的闭环（真空 + PEC），再接源；
- 选项 B：一次做到位（对偶面 J 积分 + 电荷守恒检验）。
- 事实：`EvolveE` 现在直接减 `μ0 j`；`EvolveECTRho` 在 E 更新后重算。

### D5　TM 基线

- 事实：全共形 E 会改变 2D TM 的 `r=−0.9993`（预期更接近 −1）。
- 选项 A：把旧值标为"改动前基线"，新值单独记录；
- 选项 B：以旧值为回归目标（会与更精确的格式冲突）。

### D6　AMR 的时序

- 事实：几何只在 `maxLevel()` 初始化/扩展（`WarpXInitData.cpp:1607-1633`）；`ComputeFaceExtensions`
  只对 `maxLevel()`（`WarpXFaceExtensions.cpp:557-581`）。
- 选项 A：M2a 单层，多层留 M4（与 §8.6 一致）；
- 选项 B：M2a 就做每层几何（需要把初始化/扩展改成 per-level，并处理 coarse–fine）。

---

## 8. 证据索引

| 主题 | 位置 |
|---|---|
| ECT 分配（每层） | `Source/WarpX.cpp:2683-2750` |
| ECT 初始化（仅 maxLevel） | `Source/Initialization/WarpXInitData.cpp:1607-1633` |
| 扩展链 | `Source/EmbeddedBoundary/WarpXFaceExtensions.cpp:548-633` |
| one-way / eight-ways 的 2D 限制 | 同文件 `:654-656`、`:814-816` |
| BCK 的 2D 限制 | 同文件 `:597-604`、`:613-620` |
| `ComputeSStab` | 同文件 `:84-121` |
| `GetNeigh/SetNeigh/GetNeighPtr` | 同文件 `:272-408` |
| B 三分支与 Venl | `Source/FieldSolver/FiniteDifferenceSolver/EvolveB.cpp:285-381` |
| B 循环范围 | 同文件 `:261` |
| idim=0 借用索引 | 同文件 `:299-302`、`:331-334` |
| Rho 的 3D/2D 分支 | `Source/FieldSolver/FiniteDifferenceSolver/EvolveECTRho.cpp:111-146` |
| Rho 在 E 末尾重算 | `Source/FieldSolver/WarpXPushFieldsEM.cpp:1069-1090` |
| E 更新（Yee/ECT 同一内核） | `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp:92-94`、`:171/182/193` |
| 2D 置零（边/面/flag） | `Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp:376-382`、`:465`、`:535`、`:581`、`:608` |
| 2D 掩膜（E/B） | 同文件 `:262-288`、`:323-346` |
| 2D 分量交错 | `Source/WarpX.cpp:2469-2474` |
| 2D CellSize | `Source/WarpX.cpp:3272-3286` |
| 步序 | `Source/Evolve/WarpXEvolve.cpp:646-681` |

---

## 附录 A　复现实验（摘要）

```text
/tmp/ect_test/run_te/     2D ECT, TE (Ey+Bx)，400 步，周期边界，EB 平面 z=3.9
/tmp/ect_test/run_tm/     2D ECT, TM (Ex+By)
/tmp/ect_test/run_yee/    2D Yee, TE
/tmp/ect_test/r3d_te/     3D ECT, TE，64³
/tmp/ect_test/l_te_ect/   2D ECT + 激光 (0,1,0)
/tmp/ect_test/l_tm_ect/   2D ECT + 激光 (1,0,0)
/tmp/ect_test/l_te_yee/   2D Yee + 激光 (0,1,0)
```

运行命令：

```bash
unset LD_LIBRARY_PATH
<repo>/build_2d/bin/warpx.2d inputs     # 2D
<repo>/build/bin/warpx.3d inputs        # 3D
```

判据：TE/2D ECT 下 Bx 应逐位不变、Ey 单调增；TM/2D ECT、Yee/2D、ECT/3D 下两者应随质心匀速移动。

## 附录 B　术语

- **Z 边 / X 边**：沿 z / x 方向的网格边。2D 下 Bx 在 z 边上 `(i, j+1/2)`，Bz 在 x 边上 `(i+1/2, j)`。
- **Rho**：`ECTRhofield`，面环量 $\oint E\cdot dl/S$，**不是**电荷密度（`WarpX.cpp:2738-2743` 注释）。
- **Venl**：面环量 × 面积 的累加器，用于 extended/intruded 面的记账。
- **对偶面 / 对偶边**：E 侧的几何量，位于主网格的格心/边心（半格偏移），不是 `face_areas`/`edge_lengths`。

> 本文只记录事实与待决点；若需要，我可以再单独写"每个待决点的影响面与验证方法"，或把 §4/§5
> 转成可执行的实现清单。
</content>
