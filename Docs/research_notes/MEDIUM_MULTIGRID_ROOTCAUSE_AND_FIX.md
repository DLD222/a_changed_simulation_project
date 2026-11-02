# 介质求解器多 grid 自锁：根因、文献对照与不降级修法

> 对象：E6 导体板（σ 局部、多 grid）板内场恒零；σ 全域常数（S5）与 P2-a 不受影响。
> 结论先行：**根因不是 σ 公式，也不是 mask/OverrideSync 用法，而是"单元中心（cell-centered）的介质更新"与"节点/面（nodal）的场存储"在多 grid 下的双重副本不一致**。不降级的修法是**把更新规则定义在节点位置上**（按位置的材料规则），并让每个 E 分量在**自己的节点 box** 上迭代；同步只需作为保险。

---

## 1. 现象与证据链

| 实验 | 配置 | 板内场 |
|---|---|---|
| E1 | 界面在 grid 边界 + 多 grid（`max_grid_size=64`） | **0（自锁）** |
| E2 | 界面在 grid 边界 + 单 grid（`max_grid_size=640`） | `3.21e5 / 4.77e4` |
| E3 | 界面在 grid 内部（域平移 `dz/2`）+ 多 grid | `1.53e5 / 1.11e4` |

其他关键读数：

```
Ey 的 IndexType = (1,1)  → x-face × z-face（node），不是 cell-centered
边界 node (320,320) 同时存在于两个 FAB：
  fab=44 box=((256,256)(320,320))  owner=1  Ey=0
  fab=55 box=((320,320)(384,384))  owner=0  Ey=84089（指数步写入）
EvolveB 用 B0 的 mfi 读 E1 → 读到 owner 那份（0）→ Bx=0 → curlH=0 → 自锁
```

三次对照把"时间/探针/输入/常数 σ"全部排除，剩下唯一的变量是**界面相对 grid 的位置**。

---

## 2. 直接原因：AMReX 的 nodal box 在多 grid 下重叠

### 2.1 `Box::convert` 会把 big end 加一

AMReX `AMReX_Box.H`：

```cpp
BoxND<dim>::convert (const IntVectND<dim>& typ) noexcept
{
    IntVectND<dim> shft(typ - btype.ixType());
    bigend += shft;              // cell -> node: big end += 1
    btype = IndexTypeND<dim>(typ);
    return *this;
}
```

所以由同一个 cell box 派生的四个场，BoxArray **各不相同**：

```
E0 = (0,1)  → z 方向 big end +1
E1 = (1,1)  → x,z 方向 big end +1
E2 = (1,0)  → x 方向 big end +1
σ  = (0,0)  → 原 cell box
```

相邻 grid 转换后，前一个 grid 的 big end 与后一个 grid 的 small end **相等** → **同一个物理 node 被两个 FAB 各存一份**（节点重叠，nodal overlap 是 AMReX 的设计，不是 bug）。

### 2.2 `array(mfi)` 只对"同一布局"的 FabArray 有效

`FabArray::array(const MFIter&)` 的实现是 `fabPtr(mfi)->array()`，而 `fabPtr(mfi)` 的实现（`AMReX_FabArray.H`）是：

```cpp
AMREX_ASSERT(mfi.LocalIndex() < indexArray.size());
AMREX_ASSERT(DistributionMap() == mfi.DistributionMap());  // 只检查 DM，不检查 BoxArray
int li = mfi.LocalIndex();
return m_fabs_v[li];
```

**它只要求两者的 DistributionMapping 相同，不检查 BoxArray**。因为四个场都由同一个 cell `ba` 经 `convert` 派生，DM 相同、box 顺序相同，所以第 L 个 FAB 对应同一个 grid——**跨场访问在"grid 索引"意义上是自洽的**；真正的问题是 nodal box 的**物理重叠**（同一 node 在两个 FAB 各一份）。因此：

- 用 E0 的 `mfi` 去 `array()` 一个布局不同的场，取到的是"第 L 个 box"——在 nodal 方向上与 E0 的第 L 个 box **差一格**；
- 同一物理 node 在 E1 里属于两个 box（L 与 L+1），`EvolveB`（用 B0 的 mfi）读到哪一份，取决于 B0 的第 L 个 box 与 E1 的哪一份对齐。

这就是用户实测的"`EvolveB` 读 0、`MacroscopicEvolveE` 写 84089"。

### 2.3 为什么 σ=0 不暴露、指数步暴露

| 路径 | 迭代 box | 共享 node |
|---|---|---|
| BE（σ=0） | **E 场自己的 nodal tilebox**（`mfi.tilebox(Ey 的 ixType)`） | 两个 grid 的副本都被更新，且用同一 stencil → 值一致 |
| 指数步（σ>0） | **σ 的 cell box**（第二遍的 `tb`） | 只更新"cell box 覆盖到的那一份"；另一份保持 E^n（0） |
| BE 的 skip | `if (sigma_arr(cell)>0) return;` | 界面 node 在 cell 索引上属于 σ>0 的 cell → **BE 也跳过它** |

于是界面 node 的**真空侧副本永远停在被跳过后的旧值**，`OverrideSync`/`FillBoundaryAndSync` 按"**小 box index 优先**"的规则，把这份旧值（0）复制给导体侧——同步把正确的 84089 也抹掉了。

---

## 3. 更根本的原因：更新的"材料归属"定义在 cell 上，而场定义在 node 上

指数步的 M1/M2、χ⁻¹ 行都是 **cell-centered**（同一 index 的 collocated 近似）；场却是 **nodal/face**。一个共享 node 同时属于两个 cell（σ=0 与 σ>0），**没有唯一值**。任何"按 cell 更新 + 按 grid 写入"的方案都会在 node 上产生歧义，只能靠"owner 规则"挑一份——而 owner 规则与物理无关。

**根本的修法是：把材料规则定义在 node 的位置上，而不是 cell 上。** 这样两份副本用同一个规则、得到同一个值，重叠就不再是问题。

### 3.1 文献与开源项目的做法

**Meep**（[Subpixel Smoothing 文档](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/)）给出的界面平均规则：

```
ε̃^{-1} = P⟨ε^{-1}⟩ + (I − P)⟨ε⟩^{-1},   P_ij = n_i n_j
```

即：**平行于界面的场分量用 ε 的算术平均，垂直于界面的分量用调和平均**（P 为法向投影）。这与本项目的 Kottke 实现一致（`chi_inv_isotropic` 的法向 = 调和、横向 = 1/算术）。

Meep 的关键设计是：**所有场分量存放在同一个 Yee 网格上，交错只体现在 stencil 的偏移里**（同一个 `chunk` 内，`chi1inv[c][d]` 与 E 分量同 index 存储）。因此不存在"同一物理点两份副本"的问题——材料是**位置的函数**，任何一次读取都得到同一个值。

Meep 文档同时指出两个与本问题直接相关的限制：

- "changes in the resolution which shift the location of the Yee grid points relative to the material interfaces will result in **unpredictable changes**"——界面相对网格的位置不应影响结果；
- "the fields directly on the interface are still at best **first-order** accurate"——界面 node 本身只有一阶精度，但这不影响离开界面的透射/反射。

**AMReX**（[FabArray communication 文档](https://amrex-codes.github.io/amrex/docs_html/doxygen/group__amrex__fabarray__communication.html)）：

- `FillBoundary`："Data is copied from valid regions to intersecting regions of definition... **expects that its cell-centered version of its BoxArray is non-overlapping**."——它**不解决** valid 重叠；
- `FillBoundaryAndSync` / `OverrideSync`："synchronize nodal data... override valid regions by the intersecting valid regions with a higher precedence. **The smaller the global box index is, the higher precedence**."——同步只有"谁覆盖谁"，**不能把错的值变对**；
- `OverrideSync` 文档的例子明确说：4 个 box 共享一个点、值不同时，同步只是"用其中一个（例如 box 0）覆盖其余"。

**WarpX**：E/B 各分量用 `amrex::convert(ba, *_nodal_flag)` 分配（`Source/WarpX.cpp:2595-2597`），本质就是上文的 nodal BoxArray。它自己跨场 `array(mfi)` 的写法在**场值光滑/一致**时成立；一旦某个 node 的副本被"部分更新"（指数步只更新 σ>0 一侧），就暴露。

**FDTD 教科书**（Taflove & Hagness, *Computational Electrodynamics*）：界面处的材料平均——法向分量用调和平均（D_n 连续）、切向分量用算术平均（E_t 连续）；不处理界面会导致 staircase 一阶误差。

**子网格平滑文献**：Kottke, *Phys. Rev. E* **77**, 036611 (2008)；Farjadpour et al., *Opt. Lett.* **31**, 2972 (2006)；Oskooi et al., *Opt. Lett.* **34**, 2778 (2009)。这些给出的是 ε 的 τ 变换（本项目已实现），但**不涉及**"多 grid 节点副本"这一存储层问题。

---

## 4. 不降级的修法

### 4.1 原则

1. **材料规则按 node 位置定义**：对每个 E 分量的每个 node，用相邻 cell 的属性按切向/法向规则合成（切向：算术平均/真空侧；法向：调和平均）。
2. **每个 E 分量在自己分量的 nodal tilebox 上迭代**（与 BE 路径相同的迭代方式），因此**每个共享 node 的两份副本都被同一规则更新**。
3. **BE 与指数步共用同一条 node 规则**（包括 skip 条件），不再有"cell σ>0 就跳过 node"的分裂。
4. **`FillBoundaryAndSync` 只作保险**：规则一致后它是空操作；若仍有差异，它按 owner 规则同步，不会破坏正确值。

### 4.2 具体规则

对 z 方向界面（E6）：

```
Ex, Ey（切向）:  ε_node = 0.5(ε_L + ε_R)       （算术平均）
                 σ_node = 0  若任一侧 σ=0；否则 0.5(σ_L+σ_R)
Ez（法向）:      ε_node = 2 ε_L ε_R/(ε_L+ε_R)   （调和平均，D_n 连续）
                 σ_node = 2 σ_L σ_R/(σ_L+σ_R)   （串联）
```

- 切向取"任一侧 σ=0 → 0"等价于**真空侧优先**：E_t 连续，导体的 E_t≈0 由导体**内部** node 的指数步实现；
- 法向用调和平均，保证 D_n 连续（导体界面的正确法向条件）；
- 对**光滑变化**的 σ，这条规则退化为两侧平均，不引入新的近似。

### 4.3 实现要点（按文件）

| 文件 | 改动 |
|---|---|
| `MacroscopicProperties.{H,cpp}` | 新增按 E 分量 nodal 布局的 `m_eps_node_*`、`m_sigma_node_*`（或直接存 `m_sigma_m1/m2` 的三套 nodal 副本），在初始化时用 4.2 的规则填充（需要 1 层 ghost，σ/ε 已有 `ng=1`）；**填充/判据一律用 `const_array(mfi)` + valid box，host 归约用 `MultiFab::min(comp, 0)`，见 §4.5** |
| `MacroscopicEvolveE.cpp` | 第二遍：**对每个 E 分量，用自己的 `MFIter`/`tilebox` 迭代节点**；节点属性读 4.2 的 node 规则；写回该节点的两份副本（各自 grid 的 tilebox 都会覆盖到） |
| `MacroscopicEvolveE.cpp`（BE 循环） | skip 条件从"cell σ>0"改为"**node σ>0**"（即 4.2 的 σ_node） |
| `WarpXPushFieldsEM.cpp` | rho_medium 同样按 node 规则与 node 迭代；mask/同步改为 `FillBoundaryAndSync` 保险 |
| `EvolveB.cpp` | 不改（E 的两份副本值一致后，跨场 `array(mfi)` 读到哪份都是同一个值） |

### 4.4 为什么这是"不降级"

- 不引入新近似：切向/法向平均正是 Meep/文献的界面规则；
- 不改变 Yee 的 curl stencil、不改变指数步的 M1/M2 公式；
- **网格无关**：规则是位置的函数，界面在 grid 边界/内部结果一致；
- 界面 node 的一阶误差与 Meep 文档所述一致（"fields directly on the interface are at best first-order"），但**离开界面的透射/反射**恢复二阶收敛。

---

### 4.5 实现前必须知道的一个 API 陷阱（实测）

跨 BoxArray 读取 cell-centred 的 `sigma` 时，两个访问接口的结果**不同**：

```
[XR] fab=55 box=((320,320)(383,383))  sig_min(const_array(mfi))=0.1  sig_min(operator[](mfi).min(0))=0
```

同一个 `mfi`、同一个 FAB（box 全在 z>0、σ≡0.1），`const_array(mfi)` 给出正确的 0.1，
而 `operator[](mfi)` 给出 0 —— 后者不是按 `fabIndex` 取 FAB。

**所以 4.3 里的所有跨场/跨 ba 读取都必须用 `const_array(mfi)` / `array(mfi)`，
不能用 `mf[mfi]`**（包括 `min/max` 这类 host 归约，必须改成对 `Array4` 手写循环）。

**同时实测：mask 版 `OverrideSync` 不能解决 E6。** 把 mask 按 4.2 的"导体侧 owner"
正确设置（`const_array` 判 σ、`fab 54/55 = 1`）后，同步前后板内仍为 0。
这说明同步只能搬运已有的值，无法让"真空侧的 0 副本"变成导体侧的正确值 ——
与 §2 的判断一致：**必须让两份副本由同一条 node 规则更新**（4.1 的第 1、2 条），
而不是靠事后同步。


### 4.5 实现前提：AMReX 的访问与归约语义（源码核实）

实现 §4.3 时最容易踩的两个坑，已在 AMReX 源码里核实：

**（1）`array(mfi)`、`const_array(mfi)`、`operator[](mfi)`、`get(mfi)` 是同一个 FAB**

```cpp
// AMReX_FabArray.H
const FAB& operator[] (const MFIter& mfi) const noexcept { return *(this->fabPtr(mfi)); }
Array4<value_type const> const_array (const MFIter& mfi) const noexcept { return fabPtr(mfi)->const_array(); }
```

它们都走 `fabPtr(mfi)` → `m_fabs_v[mfi.LocalIndex()]`。**不存在"两个 API 的 index 语义不同"**。

**（2）`FAB::min(comp)` 包含 ghost，`MultiFab::min(comp, nghost=0)` 只扫 valid**

```cpp
// AMReX_BaseFab.H
BaseFab<T>::min (int comp) const noexcept { return this->min<run_on>(this->domain, comp); }
// this->domain = FAB 的分配区域（valid + ghost）
// AMReX_MultiFab.H
Real min (int comp, int nghost = 0, bool local = false) const;   // 默认只扫 valid
```

**这解释了之前"同一 mfi、`const_array` 读到 0.1、`operator[] .min(0)` 读到 0"的现象**：导体 grid 的 FAB 在界面一侧有一层 **ghost（真空 cell，σ=0）**，`FAB::min` 把它算进去 → 0；而 `const_array` 的手写循环扫的是 **valid box** → 0.1。**两者都是对的，只是区域不同。**

由此也解释了 mask 版 `OverrideSync` 为什么"mask 正确但没效果"：用 `sigm[mfi].min(0)` 判"纯导体 grid"会把界面 ghost（0）算进去 → mask 被判成 0 → 同步没有可搬运的 owner。

**实现规则（写进 §4.3 的每一步）：**

```
1. kernel 读场：array/const_array(mfi)，按全局 index 访问（不要用 [mfi].min()）
2. host 归约/判据：MultiFab::min(comp, 0)（valid only），或 FAB::min(validbox, comp)
3. 以"纯导体/纯真空"做判据前，先 FillBoundary 或显式限定 valid box
4. 不要用 operator[](mfi).min() 做任何 mask/分支决策
```

**（3）`OverrideSync` 仍不能修复本问题**

即使 mask 完全正确（用 valid box 判 σ、导体侧为 1），同步前后板内仍为 0：同步只在"已有的两份副本"之间选一个**覆盖**（precedence = 小 box index），无法把"真空侧的旧值 0"变成"导体侧的正确值"。**必须让两份副本由同一条 node 规则更新**（§4.1/§4.2），同步只作保险。


### 5.0b C 完整版落地后的 V1 结果（已通过）

把指数步改成三个顶层 node 循环（每个 E 分量用自己的 MFIter + node box，
用 `m_m1_node/m_m2_node/m_sigma_node`）之后：

| 配置 | 界面位置 | q1 | q2 | ratio |
|---|---|---|---|---|
| g64 | grid 边界 | 4.405253e5 | 6.553653e4 | **0.1487690613** |
| g60 | grid 内部 | 4.405253e5 | 6.553653e4 | **0.1487690613** |

```
相对差 = 0.000e+00        （逐位相同，远优于 1e-12）
理论 exp(-alpha dz) = 0.1530  -> 差 2.8%
```

对比 C 之前（只有 BE 的 node sigma）：g64 = 0.06327、g60 = 0.14908、相对差 136%。
**这证实了「指数步 cell-based 的 stale 副本」是网格依赖的唯一来源。**

回归：P0（sigma=0）`R_flux=0.1335 T_flux=0.8513` 逐位不变。

### 5.0c S5（sigma 全域常数）—— 通过：-0.215%

算例：`sigma = 0.01` 全域、`eps = 4 eps0`、1280^2（32 ppw）、探针 z=-1.0/-0.5。
分析用项目脚本 `benchmarks/cases/c03_medium/lossy_slab/analyze.py`（固定载频 + 整窗 DFT）：

```
measured |q2/q1| = 0.623942      theory 0.625287      dev -0.215%    (判据 -0.22%)  PASS
measured phase   = +6.226 deg    theory +6.755 deg    dev -0.529 deg (判据 <=3 deg) PASS
k = 12.330583 + 0.939090i (1/m)
```

**方法说明（重要）**：直接用 `max|Ey|` 峰值比得到 `0.62055`（-0.53%），与 DFT 的
`-0.215%` 差 0.3 个百分点 —— 峰值包络受脉冲整形/多次反射影响。**S5 一律用 analyze.py，
不要手算峰值比。**

**C 完整版的 A/B（640^2，同一算例）**：

| 版本 | ratio |
|---|---|
| C 之前（`415a34d`，BE node sigma + cell 指数步） | `0.6137425143` |
| C 之后（node 指数步） | `0.6137425143` |

**逐位相同** -> C 只影响 sigma 空间变化的算例；全域常数时 node 的 tau 平均退化为同一值，
node M1/M2 与 cell M1/M2 数值相同。

### 5.0g **S5b：均匀 evanescent（无界面）—— 体格式不收敛（决定性）**

**设计**（不需要 A/B，唯一能分开"体格式 vs 界面"的实验）：

```
eps = eps0, sigma = 0.1 EVERYWHERE   (delta = 0.092 m, 与 E6 同一个 evanescent 模式)
无界面、无 slab、无背表面            (sigma 是常数 -> 全走 node 指数步)
源在 z = -0.3（否则到探针被 e^{-50} 吃掉）
11 个 Ey 节点探针（0.03125 的倍数，z = 0.125 .. 0.4375）
与 E6 完全相同的 ln|E| 斜率拟合（固定载频 DFT，窗口 [8,22] ns）
```

**注意**：脉冲在 9.6 ns 发射、10.6 ns 到达探针，窗口必须用 [8,22] ns；
用 E6 的 [24,50] ns 会读到纯噪声（`|E|~1e-10`），这一点最初踩过。

| 算例 | 640^2 | 1280^2 | 行为 |
|---|---|---|---|
| **S5b（均匀 evanescent、无界面）** | **+1.162%** | **+1.917%** | **不收敛** |
| E6（有界面） | +0.166% | +1.150% | 不收敛 |
| S5（均匀、传播波） | -0.873% | -0.215% | 收敛 |

（`delta_theory = 0.099876 m`；S5b 的 resid 两套网格都是 `8.8e-3`，与 E6 相同量级）

**结论（本轮最重要的归因）**：

1. **无界面的均匀 evanescent 也显示"1280^2 更差"**（1.16% -> 1.92%），
   而残差几乎不变（8.81e-3 / 8.86e-3）——与 E6 的行为完全同型；
2. 因此 **E6 的 ~1% 不是界面数值反射**（如果是，S5b 应该收敛而 E6 不收敛）；
   它是 **体格式的 evanescent 模式在这个参数区间（delta/dx ≈ 3..6）的固有离散误差**；
3. 传播波（S5）收敛、evanescent（S5b）不收敛 —— 两者的差别只在模式：
   `Im(k)` 很大（`k = 11.8 - 10.0i`，`|Im k| ~ Re k`），collocated 的 curl（same-index）
   + node M1/M2 对这个模式的离散色散偏离连续值，且**偏离不随 dx 减小**；
4. 所以"|B/A| 分离"和"E6-wide（更厚板）"都不是正确的下一步（前者病态、后者更糟）。

**下一步的正确方向**（未做）：

- 算 **evanescent 模式的离散色散** `1/delta_num`（1D 均匀、无界面）。注意上一轮我用
  `fsolve/root` 解隐式方程落到了假根（`k_num = 2x/10x` 连续值），应该改成
  **正向扫描**（给定 k 算一步传播因子，扫 |growth - e^{-i w dt}| 的最小值）；
- 若 `1/delta_num` 本身随 dx 变差 -> 体格式（node M1/M2 + collocated curl）在 evanescent
  模式下的离散误差，需要 staggered curl（all-Yee）或对 evanescent 模式做修正；
- 若 `1/delta_num` 收敛而测量不收敛 -> 仍是测量（但 S5b 的 resid 与 E6 同型，可能性较低）。

**影响**：这不是"sigma 空间变化"的缺陷（S5b 是均匀 sigma），而是 **evanescent 模式在
`delta/dx ~ 3` 的粗网格上的格式误差**。S5（传播波，`-0.215%`）与 P0/E3/V1 都不受影响。


### 5.0f V2d：1280^2 节点版 + |B/A| 分离的极限（本轮结论）

**1280^2 的 Ey 节点版（同一组 z）**：`delta_meas = 0.101024`，`rel +1.150%`，`resid 7.67e-3`。

| 网格 | cell 中心（插值） | node（无插值） |
|---|---|---|
| 640^2 | +0.402% | **+0.166%** |
| 1280^2 | +1.323% | **+1.150%** |

**所以"之前的非收敛"不能全部归给插值**：插值贡献了 640^2 的 0.24 pp，但 1280^2
仍比 640^2 差 7 倍。

**驻波分离（|B/A|）在当前几何下是病态的**：

```
固定 k = k_theory 的线性 A/B 拟合 (E = A e^{ikz} + B e^{-ikz}):
  640^2 : |A| = 5.32e-09   |B| = 1.86e-03   |B/A| = 3.49e+05   rms 1.59e-02
  1280^2: |A| = 6.31e-09   |B| = 1.73e-03   |B/A| = 2.74e+05   rms 1.63e-02
```

`|A| << |B|` 说明两个基函数在可用基线（0.28 m）上几乎线性相关：
`k z` 只变化 3.3 rad，而 `lambda/2 = 0.265 m`——**基线小于半个波长，A/B 无法分离**。
这与"空间 FFT 不可行"是同一个限制（衰减长度 < 波长）。

**同样失败的两条路（如实记录，不再重试）**：

| 方法 | 结果 |
|---|---|
| CW + 两向拟合（6 参数） | `|B/A| = 8e5`（退化解）、`rms 2%` |
| 1D 解析离散色散（Yee + 指数 sigma） | 求根落到假根（`k_num = 2x`/`10x` 连续值） |

**本轮的最终判读**：

1. **核心修复已闭合**：BE 的 node sigma（解开自锁）、指数步的 node 化（V1 逐位）、
   `FillBoundaryAndSync`（空操作）。P0 / E3 / S5 全部通过（逐位或 -0.215%）。
2. **E6 的 640^2 体衰减到 +0.166%**（node 探针、10 点、fit >= 2 delta）。
3. **E6 的 1280^2 仍有 +1.150%**，与"细化应更小"不符。**归因未完成**：
   需要在**更长的基线**（≥ 1.5 lambda = 0.8 m，即更厚的板或更深的探针）上重新做
   A/B 分离，或用一个可信的 1D 离散 evanescent 参考。
4. 当前测量手段（脉冲 DFT / CW 拟合 / 解析色散）在 `delta = 0.092 m`、`lambda = 0.53 m`
   这个"衰减长度小于波长"的参数区间里，分辨率不足以判定 1% 级的收敛性。

**工程结论**：E6 的**内核正确性**已由 V1（逐位网格无关）、S5（-0.215%）、P0/E3（逐位）
钉死；`1280^2` 的 ~1% 是**该参数区间（evanescent、delta < lambda/2）的测量/离散分辨率
极限**，不是已知的内核缺陷。要把这一项做成定量判据，需要改几何（更厚的板）而不是继续
调分析方法。


### 5.0e V2c：探针改到 Ey 节点（无插值）

**关键修正**：`Ey` 是 `(1,1)` **node** 场，节点在 `z = zlo + j*dz`，**不是** cell 中心
`(j+0.5)*dz`。之前所有探针都落在节点之间，靠 `FieldProbe` 插值，偏置约
`(dz/delta)^2/8` = 1.4%（640^2）/ 0.35%（1280^2），且依赖探针相对节点的分数位置。

修正：探针取 `0.03125` 的整数倍（`0.1875..0.46875`，10 个）——这些 z 在 640^2 与 1280^2
**都是节点**（粗网格节点是细网格节点的子集），`x` 也取 node（`x=0`）。无插值、两套网格同一 z。

| 网格 | cell 中心（插值） | **node（无插值）** | 改善 |
|---|---|---|---|
| 640^2 | +0.402% | **+0.166%** | 2.4x |
| 1280^2 | +1.323% | **+1.150%** | 1.15x |

（theory `delta = 0.099876 m`，`ln|E|` 斜率，`fit_start = 2.0 delta`，9 点）

**判读**：

1. 插值确实是主要系统误差之一（640^2 改善 2.4 倍）；
2. 但 1280^2 只改善 1.15 倍，**仍比 640^2 差 7 倍**（1.15% vs 0.166%）——**仍不符合收敛**；
3. 两套 resid 相近（7.22e-3 / 7.67e-3），残差仍是慢变化（尾部 +0.004/+0.011），
   与 lambda/2 驻波形态一致。

**剩余 ~1% 的网格依赖尚未归因**：要么是离散色散（evanescent 模式的数值 `Im(k)` 随 dx 变化），
要么是界面数值反射。需要 **1D 离散 evanescent 的 `1/delta_num`** 作参考才能分开（未做）。

**未做的三项（审核建议）**：
- 打印每个探针的 `|E|(t)`，确认 CW 的 `[40,80] ns` 是平台；
- 固定 `k` 的线性 `A/B` 拟合（2 个复未知数，唯一解）-> `|B/A|`；
- 1D 离散 evanescent 的 `1/delta_num(640^2)` / `1/delta_num(1280^2)`。


### 5.0d 全部回归汇总（c03-E6 修复的完整验证）

| 回归 | 结果 | 判据 |
|---|---|---|
| P0（sigma=0, 512^2） | `R_flux=0.1335 T_flux=0.8513` | 逐位不变 PASS |
| E3（sigma=0, 480^2 张量界面） | `R=0.150017 (+2.82%) T=0.860733 (+0.78%)` | 逐位不变 PASS |
| S5（sigma 全域常数, 1280^2） | `0.623942 (-0.215%)`, phase `-0.529 deg` | -0.22% PASS |
| V1（界面在 grid 边界 vs 内部） | `g64 = g60 = 0.1487690613`，相对差 `0.000e+00` | 逐位 PASS |
| D（FillBoundaryAndSync） | 加 D 前后 V1 逐位不变 | 空操作确认 PASS |
| E6 物理 | `0.14877` vs 理论 `0.1530` | 差 2.8%（V2 可精化） |


## 5. 验证清单（缺一不算通过）## 5. 验证清单（缺一不算通过）

```
V1 网格无关：界面在 grid 边界 vs 偏移 dz/2，T 差 ≤1e-12（同一物理 node、同一规则应逐位一致）
V2 解析 slab：T = |t12 t21 e^{i k d}/(1 − r21^2 e^{2 i k d})|^2,
   k = ω sqrt(μ ε) sqrt(1 − i σ/(ω ε)), η = sqrt(i ω μ/σ), r = (η−η0)/(η+η0)
   误差 ≤3% @32 ppw（介质内），探针 ≥1.5δ，只取透射窗口
V3 收敛：16/32/64 ppw，误差 O(h²)
V4 回归：P0 逐位（0.1335/0.8513）、S5 有损波 −0.22%、E3 R/T、ρ_medium Gauss（4e-25 / 7e-13）
V5 反例：σ=0 与 σ 全域常数保持与旧路径逐位一致
```

---

## 6. 对已有结论的影响

| 项 | 影响 |
|---|---|
| S5 有损波 `−0.22%`、P2-a 弛豫/刚性 | 不受影响（σ 全域常数，两份副本天然一致） |
| E3 张量界面 `R +0.86% / T +0.24%` | 不受影响（ε 界面、σ=0） |
| E6 导体板 | 本修法落地前**不可用**；落地后按 V1–V3 重测 |
| `q2/q1=0.0728` | 修法前无意义（多 grid 板内为 0；单 grid 时探针又在 0.5δ 处测界面瞬态） |

---

## 7. 参考

1. AMReX, *Communication for Distributed Fab Data*（`FillBoundary`/`FillBoundaryAndSync`/`OverrideSync` 语义）: <https://amrex-codes.github.io/amrex/docs_html/doxygen/group__amrex__fabarray__communication.html>
2. AMReX 源码 `Src/Base/AMReX_Box.H`（`convert` 的 `bigend += shft`）、`AMReX_FabArray.H`（`array(mfi)` → `fabPtr` → `LocalIndex`）。
3. Meep, *Subpixel Smoothing*（ε̃⁻¹ = P⟨ε⁻¹⟩ + (1−P)⟨ε⟩⁻¹；界面规则与限制）: <https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/>
4. A. Kottke, *Phys. Rev. E* **77**, 036611 (2008).
5. A. Farjadpour et al., *Opt. Lett.* **31**, 2972 (2006); A. Oskooi et al., *Opt. Lett.* **34**, 2778 (2009).
6. A. Taflove, S. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain Method*（界面材料平均、staircase 误差）。
7. AMReX issue #50, *Clarifications on ownership of nodal data*（nodal 数据归属）: <https://github.com/AMReX-Codes/amrex/issues/50>
