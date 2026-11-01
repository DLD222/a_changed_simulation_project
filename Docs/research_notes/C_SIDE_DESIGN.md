# C 侧设计：E 场对偶共形 Ampere（`EvolveECartesianECT`）

> 日期：2026-10-31 · 状态：**设计定稿，未实现**
> 前置：`ECT_2D_TE_REPORT.md`（B 侧已完成，L3 已通过）、`M2A_SOURCE_REVIEW_FINDINGS.md` §5（源码核查）、
> `cutcell-research.md` §8.2(b)（合并设计）
> 本文给出可直接照做的实现清单：数据结构、2D 逐分量公式、退化判据、稳定化、AMR 与 2D/3D 对照方案。

---

## 1. 为什么需要它（现状缺口）

| 侧 | 现状 | 代码事实 |
|---|---|---|
| B（Faraday） | 已是真共形（ECT 面环量） | `EvolveBCartesianECT` + `ECTRhofield` + `face_areas` |
| E（Ampere） | **仍是阶梯掩膜**：`m_eb_update_E` 只决定"更不更新"，不改变差分权重 | `EvolveE.cpp` 中 `Yee` 与 `ECT` 走同一个 `EvolveECartesian<CartesianYeeAlgorithm>` |

后果：`L2b` 中 ECT 的等效壁面之所以正确，是因为**掩膜恰好只冻结壁面节点**；一旦 surface 与网格错位
（斜切、圆弧），E 侧就退化为阶梯近似，切割单元的 Ampere 不再守恒，斜入射与绕射精度受限于 O(dx)。

C 侧目标：把 Ampere 也写成"对偶面环量"，与 B 侧对称。

---

## 2. 目标与判据

**总目标**：在 2D XZ 与 3D 上，E 更新使用对偶几何（对偶面积 + 对偶边长），PEC 上保持
`E_t = 0`，cut cell 上满足离散通量守恒；规则域必须**逐位**退回 Yee。

| 层级 | 判据 | 期望 |
|---|---|---|
| C-0 | 代码开关关闭时，结果与当前分支**逐位**相同 | `max\|d\| = 0` |
| C-1 | 规则域（EB 面移出域外）+ 开关打开：与 Yee 逐位一致 | `max\|d\| = 0`（复用 L2a 算例） |
| C-2 | 网格对齐 PEC：等效壁面位置不劣于当前（21.700 ns，误差 0.018 ns） | 反射峰误差 ≤ 0.05 ns |
| C-3 | 30° 斜面：沿束平均 `\|r+1\|` 不劣于当前 0.0124，且单点起伏收窄 | 目标 ≤ 0.01 |
| C-4 | 3D：上游 `rotated_cube` 不变（C 侧默认不改变 3D 掩膜行为时） | ≤ 1e-13 |
| C-5 | 能量守恒：真空腔 100 周期总能量漂移 ≤ 1e-10 | 与 B 侧同量级 |

---

## 3. 几何量（新增）

### 3.1 2D XZ 逐分量语义（源码核查 §5.2 已确认）

| E 分量 | 位置 | Ampere 微分 | 2D 下被切割的几何量 | 对偶面积 |
|---|---|---|---|---|
| `Ex` | x 边 `(i+½, j)` | `∂t Ex = −c²∂z By` | **对偶边长**（沿 z 的开放段） | `A_dual = 1`（每单位 y） |
| `Ez` | z 边 `(i, j+½)` | `∂t Ez = +c²∂x By` | **对偶边长**（沿 x 的开放段） | `A_dual = 1` |
| `Ey` | 节点 `(i, j)` | `∂t Ey = c²(∂z Bx − ∂x Bz)` | **对偶面积** | `A_dual = Σ(¼)·A_open` |

> 关键：`Ey` 的 `A_dual` 是**四象限各 ¼ 的开放面积之和**，不是 4 个面积直接相加（否则大 4 倍）。

### 3.2 新增数组

```
dual_face_area   : MultiFab, (1 comp), IntVect(1,1)   仅 2D 用；3D 为三维对偶面
dual_edge_length : 3 × MultiFab, IntVect = E_idim 的交错 + 半格偏移
```

2D 具体定义（`i,j` 为节点索引，`c(i,j)` 为 cell 索引）：

```
A_dual(i,j) = ¼·dx·dz·[ volfrac(i−1,j−1) + volfrac(i,j−1) + volfrac(i−1,j) + volfrac(i,j) ]
```

对偶边长（以 `Bz` 沿 x、`Bx` 沿 z 为例，采样在对偶 cell 的四个边中点）：

```
Lx(i,j) = area_frac[0](i,j)·dz      // x 法向面的开放长度，位于 (i+½, j)
Lz(i,j) = area_frac[1](i,j)·dx      // z 法向面的开放长度，位于 (i, j+½)
```

注意 AMReX 的 `area_frac[idir](i,j)` 是**该 cell 在 idir 方向的那个面**的开放度，索引约定需与
`face_areas` 一致（B 侧已经验证过这一约定，见 `ECT_2D_TE_REPORT.md` §4.1 第 5 项）。

### 3.3 退化检查（实现前先算，不写内核）

```text
规则域：volfrac ≡ 1, area_frac ≡ 1
  A_dual = dx·dz                     （与 Yee 的分母一致）
  Ex:  Lz(j±½) = dz,  A_dual = 1  →  −c²(By(j+½)−By(j−½))/dz     ✓ Yee
  Ey:  Lx = dx, Lz = dz, A_dual = dx·dz → c²[(Bx(j)−Bx(j−1))/dz − (Bz(i)−Bz(i−1))/dx]  ✓ Yee
```

**这就是 C-1 判据的来源**：若退化不成立，说明索引/半格偏移错了，不要继续加 cut cell。

---

## 4. 实现清单（按顺序）

| # | 位置 | 改动 | 验证 |
|---|---|---|---|
| 1 | `Source/Fields.H`（`warpx::fields::FieldType`） | 新增 `dual_face_area`、`dual_edge_length` | 编译 |
| 2 | `Source/WarpX.cpp` 的 ECT 分配段（现 `face_areas`/`area_mod` 同一处） | 按 3.2 分配（节点/边交错 + 1 层 ghost） | 内存/头信息 |
| 3 | `Source/EmbeddedBoundary/EmbeddedBoundaryInit.{H,cpp}` | 新增 `ComputeDualGeometry(...)`：由 `volfrac`/`area_frac` 填 3.2 的数组 | 打印 min/max，规则域应严格 `{dx·dz}`、`{dz}`、`{dx}` |
| 4 | `Source/Initialization/WarpXInitData.cpp` 的 `InitializeEBGridData` | ECT 分支调用 `ComputeDualGeometry` | 同上 |
| 5 | `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp` | 新增 `EvolveECartesianECT`，`m_fdtd_algo==ECT` 时走它（保留 `eb_update_E` 掩膜语义） | C-0/C-1 |
| 6 | 同上 | J 项：先做 `J = 0`（真空+PEC），源项留到 M2a-2 | C-1 |
| 7 | `Source/EmbeddedBoundary/WarpXFaceExtensions.cpp` | 对偶 `S_stab`（把 3.2 的 `A_dual` 抬到 `0.5·max{邻接长度×间距}`），BCK 兜底 | 小胞元扫描 |

**开关**：建议加 `algo.ect_conformal_E = 0/1`（默认 0），先保证 C-0（关掉时逐位不变），再打开做 C-1/C-2。

### 4.1 2D `Ey` 的离散（关键公式）

```
A_dual(i,j)·(Ey_new − Ey_old) =
    c²·dt·[ (Bx(i, j)·Lz(i,j) + Bx(i−1, j)·Lz(i−1,j)) / 2
          − (Bx(i, j−1)·Lz(i,j−1) + Bx(i−1,j−1)·Lz(i−1,j−1)) / 2
          − (Bz(i, j)·Lx(i,j) + Bz(i, j−1)·Lx(i,j−1)) / 2
          + (Bz(i−1, j)·Lx(i−1,j) + Bz(i−1,j−1)·Lx(i−1,j−1)) / 2 ]
```

（上游/下、左/右各两段，`/2` 来自每段只占 cell 面的一半；规则域化简后与中心差分一致。）

> 符号与左右归属必须在规则域上用"逐位一致"钉死；建议先写一个 `A_dual ≡ dx·dz`、
> `L ≡ 常数` 的临时版本，跑 L2a，确认 0 差异后再接入真实几何量。

---

## 5. 稳定化（对偶侧）

与 B 侧对称，但对象不同：

| 侧 | 小量 | 稳定化 |
|---|---|---|
| B | 小 `face_areas` | ECT 扩展/借用（已有）+ BCK |
| E | 小 `A_dual`（Ey）、小 `L`（Ex/Ez） | 先做 BCK：`A_dual < A_stab` 时抬到 `A_stab`；`L` 同理 |

`A_stab = 0.5·max{ 相邻两段的 L·间距 }`，与 `ComputeSStab` 同构。
**不建议**直接搬 B 侧的 `Venl` 记账：E 侧没有"环量累加器"的等价物，先 BCK-only 跑稳定性扫描，
再决定是否引入借用。

---

## 6. AMR 多层几何（当前实现的硬限制）

代码事实（`ECT_2D_TE_REPORT.md` §4.3、`M2A_SOURCE_REVIEW_FINDINGS.md` §1.1）：

- `ComputeEdgeLengths` / `ComputeFaceAreas` / `ComputeFaceExtensions` / `MarkExtensionCells`
  只在 `lev == maxLevel()` 调用（`InitializeEBGridData` 内）；
- `face_areas` 等数组虽然按层分配，但**低层从未初始化**；
- `warp` 侧 `EvolveB` 用 `face_areas[lev]`，在 `lev < maxLevel` 时读到的是未初始化/零值。

要做多层，需要：

1. 把 `InitializeEBGridData` 的 ECT 分支改成 **每层** 调用（`lev == maxLevel` 的条件去掉），
   并按该层的 `cell_size`/`eb_fact` 计算；
2. `ComputeFaceExtensions` 的 `maxLevel()` 限制（`WarpXFaceExtensions.cpp`）改成 per-level；
3. 粗-细界面：`area_mod`/`Venl` 在 `lev` 与 `lev+1` 的重叠区需要一致化（插值或重新计算 + `FillBoundary`）；
4. `MarkUpdateBCellsECT` 的 B 掩膜同样按层；
5. 验证：两层 AMR + 斜 PEC，比较"细层覆盖区"与"单层细网格"解一致（≤1e-8），再测守恒。

**建议顺序**：C 侧单层 → L3 复测 → AMR 每层几何 → AMR 多层求解。

---

## 7. 2D/3D 逐点一致性对照（决策 #11）

设置（硬约束）：

```
3D 域：z 与 x 同 2D；y 方向 1 层网格 + 周期边界；dy = 1（与 per-unit-y 约定对齐）
同 dx、dz、dt、同初值表达式（把 2D 的 parser 公式原样用在 3D，y 不出现）
同诊断：探针 (x, z) 相同、y = 0
容差：逐点相对差 ≤ 1%（先看量级；若格式相同应远小于此）
```

判据脚本：读两次运行的 `fp1.txt`，逐列比较 `Ey/Bx`（TE）与 `Ex/By`（TM），输出 `max|d|/peak`。

风险：3D 的 ECT 走的是 3D 面环量（`Rhox/Rhoy/Rhoz` 全量），2D 走 per-unit-y 约化；
两者在 **cut cell** 上是否逐点一致**没有先验保证**，只在规则域必然一致（B 侧已验证）。
因此该对照的预期是"规则域一致、cut cell 下差异 O(dx)"，需要在文档里明确，不要把 1% 当硬门槛。

---

## 8. 工作量与风险

| 项 | 估计 | 风险 |
|---|---|---|
| 对偶几何（2D） | ~150 行 + 自检 | 索引/半格偏移错 → 退化不平；用 L2a 逐位判据兜底 |
| `EvolveECartesianECT`（2D Ey） | ~120 行 | 符号/左右归属；同上 |
| Ex/Ez 对偶边长（2D TM） | ~100 行 | TM 基线 `r=−0.9994` 会变，需重记基线 |
| 稳定化 | ~80 行 | 小胞元扫描缺失会漏稳定性问题 |
| 3D | ~250 行 | 三个分量对称性易错，建议先 2D |
| AMR per-level | ~120 行 | 粗-细一致性 |

**风险最高的不是公式，而是索引约定**（`area_frac` 的左/右、`IntVect` 半格、`k` 维度在 2D 的
哑元）。B 侧的经验是：先用规则域逐位判据把索引钉死，再加 cut cell。

---

## 9. 下一步（可直接执行）

1. 写 `ComputeDualGeometry`，只打印规则域的 min/max，确认 `A_dual = dx·dz`、`L = dz/dx`；
2. 写死几何（`A_dual`、`L` 常数）的 `EvolveECartesianECT`，跑 L2a，要求 `max|d| = 0`；
3. 接入真实几何量，重跑 L2a（仍须 0）与 L2b（等效壁面不劣于 21.700 ns）；
4. 跑 L3，比较沿束平均 `|r+1|`（当前 0.0124）；
5. 打开 BCK 对偶稳定化，做小胞元扫描；
6. 最后做 AMR 每层几何与 2D/3D 对照。
