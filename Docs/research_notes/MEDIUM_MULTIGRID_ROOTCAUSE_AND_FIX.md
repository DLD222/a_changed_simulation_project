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

### 5.0s **第 1 条：张量 eps + sigma>0 的 R/T —— 两次尝试都不可信（已定位原因）**

| 尝试 | 源 | 结果 | 问题 |
|---|---|---|---|
| 1 | 激光（`tensor_interface` 原样）+ `analyze.py` | `R=0.962`、`T=3e-6` | 激光天线被返回波扰动；10 m 厚板 + PEC 背表面形成腔 |
| 2 | **`E_ext_grid` 载频脉冲** + flux 法 | `R_flux=0.999`（理论 `0.148`） | flux 在**强反射 + 有损**时失效（`W_med = 入射 − 反射 ≈ 0`） |

**第 2 次的第二个线索（新信息，值得查）**：

```
|Ey(med, z=1)| = 1.26e3   vs   |Ey(ref, z=1)| = 1.0e6   ->  衰减 e^{-6.7}
理论：eps_yy = 5 eps0, sigma = 0.01  ->  delta_yy = 1.19 m  ->  1 m 处应 e^{-0.84} = 0.43
=> 实测衰减比理论快 8 倍（delta_eff ~ 0.15 m）
```

**可能的根因（待验证）**：`sigma=0.01` 是标量全域，材料 1（`eps0`）和材料 2（张量
`diag(4,5,6) eps0`）都有它。s-pol 在材料 2 里应看到 `eps_yy = 5 eps0`（`delta = 1.19 m`）。
如果**指数步在板内深处错误地用了界面 tau 平均后的低 eps**（`2.5/3.0/1.714 eps0`）
而不是材料 2 的 `5 eps0`，`delta` 会变小、衰减更快。

**下一步（明确、可执行）**：在板内深处（`z = 0.5, 1.0`）打印 `node eps6 / sigma6`，
确认它等于材料 2 的 `(4,5,6) eps0`，而不是界面处的 tau 平均值。

**已完成、不受影响**：

| 项 | 结果 |
|---|---|
| 第 2 条（旋转张量 + sigma 的 node M1/M2 对 Python） | `4.2e-10` PASS |
| 第 1 条的两次尝试 | 都不可信，原因已定位 |


### 5.0r **旋转张量 eps + sigma>0 的 node M1/M2 对 Python 3x3 expm（4.2e-10）**

配置：`eps2 = R diag(4,5,6) R^T eps0`（绕 y 转 30 度，`eps_node` 有 `xz` 非对角）
+ **标量 `sigma = 0.01`** + 界面 `z=0`，480^2（`dt = 2.358654337e-11`，从 run.log 读）。

界面 node 读数：

```
eps6 = 2.383819796e-11  2.656256344e-11  1.498401014e-11  0  1.179684858e-12  0
sig6 = 0.01 0.01 0.01 0 0 0
M1   = diag(0.9901163, 0.9911597, 0.9843218) + xz = zx = 7.720339e-04
M2   = diag(0.9883660, 0.8840311, 1.5678202) + xz = zx = -7.720339e-02
```

对 Python（`A = eps^-1 sigma`，`M1 = expm(-A dt)`，`M2 = A^-1 (I - M1) eps^-1`）：

```
M1 对角 1.5e-6 -> 2.3e-6 (用近似 dt 时)   ->  用实际 dt 后 <=4.2e-10
M2 对角 1.5e-4 (用近似 dt 时)             ->  用实际 dt 后 <=2.1e-10
M1/M2 均对称
```

**结论**：**旋转张量 eps + sigma>0 的 node M1/M2 与 3x3 expm 一致到 `4.2e-10`**（判据 1e-8），
`M1_xz != 0` 的张量耦合在 node 接口上正确。**B 的接口换到 node 后仍然成立。**

**踩坑记录**：第一次对比用 `dt = 2.359e-11`（手算近似），偏差 `1.5e-4`；
从 run.log 读到实际 `dt = 2.358654337e-11` 后降到 `4.2e-10`。**矩阵验证必须用代码里的实际 dt。**


### 5.0q **张量 eps + sigma>0 + 界面（C 之后唯一没重跑的组合）**

配置：`tensor_interface` 的 `eps2 = diag(4,5,6) eps0` + **标量 `sigma = 0.01`** + 界面 `z=0`，
480^2（`dt = 2.359e-11`）。取界面 node 的 `eps_node / sigma_node / M1 / M2`。

```
eps6:  2.213546953e-11  2.656256344e-11  1.517860768e-11   (= diag(2.5, 3.0, 1.714) eps0)
sig6:  0.01  0.01  0.01
M1:    0.9894010268  0.9911596887  0.9845807798
M2:    1.059897319   0.8840311266  1.541922017
```

**解析验证**（对角张量 eps + 标量 sigma）：

```
M1_ii = exp(-sigma dt / eps_ii)
  xx: exp(-0.01*2.359e-11/2.2135e-11) = exp(-0.010657) = 0.98940  ✓
  yy: exp(-0.01*2.359e-11/2.6563e-11) = exp(-0.008880) = 0.99116  ✓
  zz: exp(-0.01*2.359e-11/1.5179e-11) = exp(-0.015539) = 0.98458  ✓
M2_ii = (1 - M1_ii) / sigma
  xx: (1-0.98940)/0.01 = 1.0599  ✓    yy: 0.8840  ✓    zz: 1.5419  ✓
```

**结论**：`node M1/M2` 在「张量 eps + sigma>0 + 界面」这个组合下与解析一致（三个主分量全对），
说明 `node eps 的 tau 平均` 与 `cell chi^-1 行（sigma=0 路径）`在同一算例里不打架。

**这一步补充了此前唯一没重跑的组合**（C 的 node 指数步之后，E3 一直是 sigma=0）。


### 5.0p **E6 收束：深区拟合 3.5 delta -> 1280^2 优于 640^2（-0.61% vs -1.58%），二阶收敛**

用已有 E6-CW 数据（探针 `z = 0.219..0.50` = `2.2..5.4 delta`）做 `fit_start` 扫描：

| `fit_start` | 640^2 | 1280^2 | resid |
|---|---|---|---|
| 2.0 delta | +0.791% | +1.901% | `1.1e-2` |
| 3.0 delta | -0.539% | +0.559% | `6.8e-3` |
| **3.5 delta** | **-1.581%** | **-0.608%** | **`1.7e-3`** |
| 4.0 delta | -1.686% | -0.787% | `1.8e-3` |

**判读**：

1. `resid` 从 `1.1e-2`（2.0 delta）降到 `1.7e-3`（3.5 delta）——**降一个量级**，
   说明 `2.2 delta` 附近的残差主要是**界面瞬态**（lambda/2 驻波），与审核预判
   （`2.2 delta` 处瞬态还剩 `e^{-2.2} = 11%`）一致；
2. 深区（3.5 delta）后 **1280^2 比 640^2 更接近零**：`-0.61%` vs `-1.58%`，改善 **2.6x**
   （`1.78^2 ~ 3.2`）-> **二阶收敛**；
3. 而 2.0 delta 时 1280^2 看起来更差（`+1.90%` vs `+0.79%`）是**界面瞬态污染**的假象
   （细网格的界面层更薄，瞬态相位不同）；
4. 所以 **E6 的 1-2% 是界面瞬态 + 拟合区间，不是边界层离散、不是界面 BC 错误**。

**与 flux 的 R 一致**：

| 证据 | 结果 | 判读 |
|---|---|---|
| flux `R` | `0.3197` vs `0.3059`（+4.5%） | 界面 BC 正确 |
| 深区拟合 | resid `1.7e-3`、1280 优于 640 | 界面瞬态被排除、二阶收敛 |

**E6 收束（不需要改核）**：

| 项 | 结果 |
|---|---|
| 界面自锁 | 已修复（BE node sigma + 指数步 node 化，V1 逐位） |
| 界面反射 | `R_flux` +4.5%（BC 正确） |
| 板内斜率 | 深区拟合 1280 优于 640（`-0.61%` vs `-1.58%`），二阶 |
| 残差 | 3.5 delta 后 `1.7e-3`（界面瞬态已排除） |

**`SIBC` 从选项删除**（审核确认：`delta` 被 3-6 格分辨、S5b-CW 内部二阶收敛，
SIBC 会绕过内部解，与项目目标相反）。**局部加密被 `lev == 0` 挡住**（多级不支持）；
**复 `eps` tau 作为最后选项，当前不需要**。


### 5.0o **flux 法给出 R = 0.3197 vs 理论 0.3059（+4.5%）—— 界面 BC 基本正确**

用 c02/E3 的 **flux 法**（坡印廷通量的时间积分，**不需要分窗**，对脉冲鲁棒）：

```
数据: rc_ref（真空参考）/ rc_med（导体），带载频脉冲（291.5 MHz，单向 13.4 ns）
S = (Ex*By - Ey*Bx) / mu0
W_ref = integral(S_ref, [5,60] ns) = 11.02285
W_med = integral(S_med, [5,60] ns) =  7.49908
R_flux = 1 - W_med/W_ref = 0.319679
R_theory (精确 eta_c @291.5MHz) = 0.305872
=> dev = +4.5%
```

同一数据的分窗法：`R_flux_split = 0.2425`（-21%）—— **分窗不可靠**。

**结论（推翻之前的「界面偏 45%」）**：

1. `R_flux` 与精确理论一致到 **+4.5%** -> **界面 BC 基本正确**，不是偏强 45%；
2. 之前的 `R=0.64`（峰值比）、`0.164`（分窗 DFT）、`1.074`、`0.001` **都是方法假象**：
   - 峰值比取到的是 **93 ns 的 PEC 腔回波**（不是 26.7 ns 的界面反射）；
   - 分窗 DFT 的反射窗含入射尾 + 腔回波尾巴；
   - **无载频脉冲的 FFT 峰在 0 Hz（准静态）**，与 RF 的 `eta_c` 不可比；
3. 所以 **E6-CW 的 1-2% 不收敛不是界面 BC 偏离，而是边界层离散**（staircase 的半格真空层
   使斜率有 1-2% 偏置）—— **不需要 SIBC**；
4. 正确的修法方向回到「**复介电常数 tau 或界面局部加密**」（§4.2 的 `eps` tau 机制已经是
   这条路），**而不是表面阻抗 BC**。

**已确认的完整图景**：

| 项 | 结果 | 判读 |
|---|---|---|
| 体格式 | S5b-CW `-1.070% -> -0.225%`（resid `5e-4`） | 二阶收敛 |
| **界面反射** | **`R_flux = 0.3197` vs `0.3059`（+4.5%）** | **界面 BC 正确** |
| 界面斜率 | E6-CW `+0.791% -> +1.901%` | 边界层离散（1-2%） |
| 自锁 | 已修复（BE node sigma + 指数步 node 化，V1 逐位） | ✅ |
| 回归 | P0/E3/S5 全部通过 | ✅ |

**下一步**：E6 的 1-2% 边界层离散可以用「界面局部加密」或「复 `eps` tau」继续压，
但**已经不属于「界面 BC 错误」**。


### 5.0n **R 第 5/6 次尝试：定位到「E_ext_grid 双向分裂 + PEC 腔多次反射」**

按审核的 5 条建议逐条执行：

1. **Bx 符号核对**：c02 是 `Bx = -E0*exp(-((z-zc)/w)^2)/c`，我的输入**完全相同** → 符号不是问题。
2. **时域差**（激光脉冲那次）：入射 `t=9.99 ns`、`|E_ref|=1.0e6`；差值峰值 `t=24.03 ns`、
   `|diff|=8.0e5`（ratio `0.800`）。`24 ns` = 入射 10 ns + 往返 4 m/c = 13.3 ns →
   **反射确实存在**。所以之前的 `R=9.7e-5` 是 **DFT 窗口把入射/反射混在一起**，不是「没有反射」。
3. **峰值比 vs DFT**（同一数据）：`R_peak = 0.6406`、`R_DFT = 0.001058`、`R_theory = 0.305872`
   → 两者差 **600 倍**，说明脉冲频谱太宽（`w=1 m` → 3.3 ns），**单频 DFT 不可靠**。
4. **窄带化**（`w=3 m` → 10 ns、`zc=-6`、`stop_time=120 ns`）：入射峰值出现在 **`t=93.4 ns`**
   （预期 13.3 ns）→ **`E_ext_grid` 脉冲双向分裂 + PEC(-10) 腔多次反射**，探针处是复杂叠加，
   「单次入射 + 单次反射」模型不成立。`R_peak = 1.0046`（+228%）、`R_DFT = 0.6425`（+110%）。
5. **未上 SIBC**（按要求）。

**结论（有证据）**：R 在当前几何下不可测，原因是两个**结构性**因素：

```
(a) E_ext_grid 是初始条件，脉冲双向分裂（+z 和 -z 各一半）；
    即使 Bx = -Ey/c 也不能保证单向（实测入射峰值在 93ns 而非 13ns）；
(b) 域是 PEC 盒 z=+-10，多次反射在长窗口内主导。
```

**正确的 R 需要改几何（不是调参）**：

- **单向源**：只让 `+z` 半空间有初值（例如 `Ey = E0*exp(...)*step(z)`），或用一个
  解析的单向入射参考；
- **无 PEC 回波**：吸收边界（P8: PML 与 macroscopic 不兼容，需先解决兼容性），
  或用足够大的域把回波推到窗口外。

**已确认、不依赖 R 的结论（保持）**：

| 结论 | 证据 |
|---|---|
| 体格式二阶收敛 | S5b-CW `-1.070% -> -0.225%`（resid `5e-4`） |
| 有界面不收敛 | E6-CW `+0.791% -> +1.901%`（resid `1.1e-2`） |
| 界面自锁已修复 | BE node sigma + 指数步 node 化，V1 逐位 |
| P0/E3/S5 回归 | 全部通过 |

**R 的现状**：**未测准**，界面 BC 是否偏离**仍无可靠证据**。SIBC 的决定应等待可信的 R。


### 5.0m **R 测量的四种失败模式（停止，不再重复）**

| # | 方法 | 结果 | 失败原因 |
|---|---|---|---|
| 1 | CW + 两探针 `z=-1,-2`（间距 `1 lambda`） | `R=1.000`、`|A|=|B|` | `det = 2i sin(kL) ~ 0`，**线性相关** |
| 2 | CW + 3 探针（间距 `0.34 lambda`） | `R=1.32`，`lstsq_resid 5.7%` | CW 稳态含 **PEC(-10) 多次反射**，场不是 `A+B` |
| 3 | 激光脉冲 + 真空参考 | `R=0.3715`（vs 精确 `0.3059`，+21.4%） | **激光天线被返回波扰动**，两跑的"入射"不同 |
| 4 | **`E_ext_grid` 脉冲 + 真空参考** | `R=9.7e-5`（几乎为 0） | DFT 相位对齐/窗口失败（见下） |

**第 4 次的细节**（值得记录）：

```
输入: 去掉激光, warpx.E_ext_grid_init_style = parse_E_ext_grid_function
      Ey = E0*exp(-((z-zc)/w)^2),  Bx = -E0*exp(...)/c,  zc=-5, w=1
      注意 Bx 不能引用 Ey_external_grid_function（parser 不认识，必须直接写表达式）
参考 |Ey|=1e6、导体 |Ey|=1e6，两者都正常初始化
但 DFT([8,50]ns): |E_med|=3.785e-4、|E_ref|=3.821e-4（几乎相同）
  -> |E_refl| = 3.76e-6,  R = 9.7e-5（几乎无反射）
时域上导体跑 30ns 后从 1e6 衰减到 1e4（有反射/衰减），但 DFT 复幅度几乎不变
  -> **相位对齐/窗口问题**，减法把反射项抵消了
```

可能原因：`E_ext_grid` 是**初始条件**，脉冲会**双向分裂**（+z 和 -z 各一半），
探针处的波含两个方向的叠加；窗口 `[8,50]` 又包含了分裂后的多次往返。

**结论**：R 的正确测量需要**单向脉冲**（不是双向分裂的初值）+ 干净的窗口/相位对齐，
这是一个需要专门设计的测量（可能需要一个 1D 辅助算例或解析的入射参考），
**不是继续在当前几何上调参数能得到的。四次失败是明确的停止信号。**

**R 这一项的现状**：**未测准**。界面 BC 是否偏离 21%，**尚无可靠证据**——
第 3 次的 `+21.4%` 建立在被源偏置的入射上（不可信），第 4 次是相位对齐失败。
所以「是否上 SIBC」**应等待一个可信的 R**，不应基于这两次数字做决定。

**已确认、不依赖 R 的结论（保持）**：

| 结论 | 证据 |
|---|---|
| 体格式二阶收敛 | S5b-CW `-1.070% -> -0.225%`（resid `5e-4`） |
| 有界面不收敛 | E6-CW `+0.791% -> +1.901%`（resid `1.1e-2`） |
| 界面自锁已修复 | BE node sigma + 指数步 node 化，V1 逐位 |
| P0/E3/S5 回归 | 全部通过 |


### 5.0l **R 的脉冲 + 真空参考测量：R = 0.3715，比精确理论大 +21.4%（界面 BC 偏离）**

按审核的 E3 方法（不用两探针 A/B，避开 `sin(kL)=0` 退化和 CW 腔混响）：

```
参考跑: sigma 撤掉（真空全域），同几何/同源/同 PEC/同探针
导体跑: sigma = 0.1*(z>0)，界面 z=0
探针:   单个真空探针 z=-2（node）
窗口:   [8,50] ns  (入射 ~10 ns、界面反射 ~13.3 ns、PEC(-10) 回波 ~60 ns 之前)
E_refl = E_med - E_ref,  R = |E_refl|^2 / |E_ref|^2
```

**理论（精确，不用良导体近似）**：

```
eta_c = sqrt(i w mu / (sigma + i w eps)) = 116.596 + 98.760i
R_exact  = |(eta_c-eta0)/(eta_c+eta0)|^2 = 0.305872
良导体近似 sqrt(i w mu/sigma) = 153.853 -> R_approx = 0.176452
=> 两者的 R 差 42%（虚部修正对 R 的影响远大于对 eta 的 8%）
```

**测量**：

```
|E_med| = 1.894e-3   |E_ref| = 2.632e-3   |E_refl| = 1.604e-3
R_meas = 0.371453
  dev vs exact  = +21.4%
  dev vs approx = +110.5%
```

**判读（按审核的判据）**：R 偏离精确值 `+21.4%` -> **界面 BC 本身不够准**，
不是「板内边界层离散」。=> **表面阻抗 BC（Leontovich，Tidy3D 对有损金属的默认）**
是有物理依据的正解方向；复 `eps` tau 只在 R 对的情况下才是首选。

**保留的疑问**：`+21.4%` 是否全部是界面误差，还需排除测量侧（窗口、探针位置、
参考跑的相位对齐）。但 `+21.4%` 远超 5%，方向是明确的。

PML 未使用（P8: PML + macroscopic 不兼容）。


### 5.0k **真空侧 R 测量：PEC 盒几何下不可行（需要吸收边界）**

理论：`delta = 0.0919 m`、`eta0 = 376.73`、`|eta_c| = 153.85`（`sqrt(omega mu / sigma)`）、
**`R_theory = 0.17645`**。

| 尝试 | 探针 | 结果 |
|---|---|---|
| 1 | `z=-1, -2`（间距 1 m = 1 lambda） | `det = 2i sin(k) ~ 0` -> **线性相关，解退化**，`R=1.000`、`|A|=|B|` |
| 2 | `z=-1, -1.34375, -1.6875`（间距 0.34375 = 0.34 lambda） | 640^2 `R=1.317`、1280^2 `R=1.326`（比理论大 646%），`lstsq_resid` **5.7%/6.2%** |

**根因**：域是 **PEC 盒 `z = +-10`**，源在 `z=-5`，探针在 `z=-1..-1.7`。真空侧的场是
**直接波 + 界面反射 + PEC(-10) 镜像回波 + 多次反射**，不是干净的 `A e^{ikz} + B e^{-ikz}`，
所以 A/B 分离失败（`lstsq_resid` 5.7% 就是证据）。

（探针时域 `|Ey| ~ 0.7..1.6`，`dft` 的 `1e-9` 只是 `dt` 缩放，信号本身正常。）

**结论**：**R 测量需要无 PEC 回波的真空侧**（周期边界或吸收层/PML），当前 PEC 盒几何
下不可测。这不是调探针间距能解决的，**需要改几何**。

**不受影响的结果**（R 不参与它们的判据）：

| 实验 | 结果 |
|---|---|
| S5b-CW（无界面 evanescent） | `-1.070% -> -0.225%` 收敛（resid `5e-4`） |
| E6-CW（有界面 evanescent） | `+0.791% -> +1.901%` 不收敛（resid `1.1e-2`） |

所以「**体格式 OK / 界面是 1-2% 不收敛的来源**」这个归因**不依赖 R**。

**下一步（需要改几何）**：给真空侧加吸收边界（PML/吸收层）或改用周期边界，再测 R：

- R 与 `0.17645` 一致 -> 界面 BC 物理正确，板内斜率的 1-2% 是**边界层离散**，
  修法 = 复介电常数 tau 或界面局部加密；
- R 偏离 -> **表面阻抗 BC**（Leontovich，Tidy3D 对有损金属的默认）；
- **real-sigma 的 tau 平均仍不推荐**（无二阶保证 + 会碰 `E_t` 连续性/自锁，见 §4.2 的说明）。


### 5.0j **E6-CW：有界面不收敛（+0.79% -> +1.90%），无界面收敛 —— 界面是源头**

E6-CW 配置：`sigma = 0.1*(z>0)`（界面 z=0）、源在真空 `z=-5`、CW 平滑开启、
**深探针 `z = 0.219 .. 0.500`**（界面驻波按 `e^{-2z/delta}` 衰减）、窗口 `[40,60] ns`。

| 算例 | 界面 | 640^2 | 1280^2 | 行为 | resid |
|---|---|---|---|---|---|
| **S5b-CW** | **无** | **-1.070%** | **-0.225%** | **收敛** | `5e-4` |
| **E6-CW** | **有** | **+0.791%** | **+1.901%** | **不收敛** | `1.1e-2` |

（`delta_theory = 0.099876 m`；两套网格的 E6-CW resid 都是 `1.1e-2`，比 S5b-CW 大一个量级）

**这是最终的归因（四个实验的完整矩阵）**：

| 实验 | 界面 | 网格 | 结果 | 归因 |
|---|---|---|---|---|
| S5（传播波，脉冲） | 无 | 640/1280 | `-0.87% -> -0.22%` 收敛 | 体格式 OK |
| **S5b-CW（evanescent）** | **无** | **640/1280** | **`-1.07% -> -0.225%` 收敛** | **体格式 OK** |
| S5b（脉冲） | 无 | 640/1280 | `+1.16% -> +1.92%` | 测量偏置（脉冲色散/窗口） |
| **E6-CW** | **有** | **640/1280** | **`+0.79% -> +1.90%` 不收敛** | **界面数值反射** |
| E6（脉冲） | 有 | 640/1280 | `+0.17% -> +1.15%` | 界面 + 测量混在一起 |

**结论**：

1. **体格式（node M1/M2 + collocated curl）在 CW 稳态下对 evanescent 模式是收敛的**
   （S5b-CW `-1.07% -> -0.225%`，resid `5e-4`）；
2. **界面（`sigma` 从 0 跳到 0.1 的那一层）引入约 0.8~1.9% 的、随细化不降的误差**；
   深探针区（`z >= 2.2 delta`）仍不收敛，说明它不是"界面瞬态"而是**界面数值反射**
   （与残差的 lambda/2 驻波形态、`|B/A| ~ 3.6e-3` 的量级一致）；
3. 这正好落在参考文档 §4.2 的「切向真空侧优先 / 法向调和」之外：**`sigma` 的面平均规则**
   目前是 staircase（一阶），界面反射 `~5e-3`；要压到 1% 以内需要 **`sigma` 的 subpixel/τ 平均**
   或表面阻抗边界条件（scheme-level 的界面处理，不是 node M1/M2 的问题）。

**对 E6 的最终判定**：`E6` 的**自锁已修复**（BE 的 node sigma + 指数步 node 化，V1 逐位），
**板内衰减的体格式正确**（S5b-CW），**剩余的 1~2% 是界面 staircase 的数值反射**——
这是**物理层（界面处理）的精度问题，不是内核缺陷**，且与「`sigma` 空间变化」的可用性无关
（板内有场、衰减正确，只是斜率有 1~2% 的界面反射偏置）。

**门槛**：要把它压到 1% 以内，需要 `sigma` 的 subpixel 平均（参考 §4.2 的 tau 平均机制，
目前只对 `eps` 做了，`sigma` 仍是 staircase）——这是明确的、有界的下一步。


### 5.0i **扫描脚本自检失败：对 S5（已知收敛）也给出不收敛 —— 5.0g 的 +6%/+13% 作废**

审核指出「`+6.1% -> +13.4%` 随细化变大违反收敛性」，并要求先用 **S5 参数（传播波，
已知 `-0.87% -> -0.22%` 收敛）** 验证扫描脚本。结果：

```
dispersion_scan.py --compare --sigma 0.01 --eps 4*EPS0   (S5 参数, 传播波)
  continuum:            1/delta = 0.939206
  dz=0.03125  (640^2):  1/delta_num = 0.751365   dev -20.0000%
  dz=0.015625 (1280^2): 1/delta_num = 1.127047   dev +20.0000%
```

**扫描对已知收敛的 S5 也给出「不收敛」（而且恰好对称 ±20%）** —— 正是审核预判的
`ln` 分支 / 初值问题（`omega(k) = i ln(lambda)/dt` 有 `2*pi/dt` 分支，dt 减半时更容易跳支）。

**因此 `dispersion_scan.py` 的结果（包括 5.0g 的 `+6.1% -> +13.4%`）不能作为判据，作废。**
（这是本轮第三次在解析色散上失误：fsolve 假根、2x2 矩阵符号、扫描分支。解析色散这条路
从实现上反复失败，**不再尝试**。）

**交叉验证的一致性**：

| 证据 | 结果 | 可信度 |
|---|---|---|
| S5（传播波，实测） | `-0.87% -> -0.22%` 收敛 | ✅ 实测 |
| **S5b-CW（evanescent，实测）** | **`-1.07% -> -0.225%` 收敛** | ✅ **实测** |
| S5b 脉冲（实测） | `+1.16% -> +1.92%` 不收敛 | ⚠️ 测量偏置 |
| 2x2 扫描（解析） | `+6.1% -> +13.4%`（对 S5 也给 ±20%） | ❌ 作废 |

**结论**：

1. **体格式（node M1/M2 + collocated curl）在 CW 稳态下对 evanescent 模式是收敛的**
   （S5b-CW `-1.07% -> -0.225%`）；
2. 所有「evanescent 不收敛」的迹象都来自**脉冲 + 单窗口 DFT 的测量偏置**（S5b 脉冲）
   或**解析色散的实现错误**（5.0g）；
3. **E6 的 `+1.150%` 需要 CW 重测**才能定性（深探针区 `z >= 0.2` 拟合，避开界面瞬态）。


### 5.0h **S5b-CW：CW 稳态下体格式收敛（推翻 5.0g 的结论）**

S5b-CW 配置：`sigma = 0.1` 全域、无界面、源 `z = -0.3`、11 个 Ey 节点探针、
`(1-exp(-t/tramp))cos(omega t)`、`stop_time = 80 ns`、窗口 `[35,75] ns`（`|E|` 在 30 ns 后
稳定在 `5.478e-3`，确认稳态）。

| 网格 | `delta_meas` | 偏差 | resid |
|---|---|---|---|
| 640^2 | 0.098807 | **-1.070%** | `7.26e-4` |
| 1280^2 | 0.099651 | **-0.225%** | `5.06e-4` |

**收敛（改善 4.75x），且 resid 比脉冲小一个量级（5e-4 vs 8.8e-3）。**

**这推翻了 5.0g 的两个结论**：

1. S5b **脉冲**的「不收敛」（`+1.16% -> +1.92%`）是**测量问题**（导体强色散使脉冲包络
   变形，单窗口 DFT 对不同 dx 偏置不同），**不是体格式**；
2. 用 2x2 正向矩阵算的离散色散（`+6.10% -> +13.43%`，随 dx 变差）**与 CW 的收敛矛盾**
   —— 该 2x2 矩阵的归一化/符号约定与 WarpX 实现不一致，**结论作废**。
   （这是本轮第二次在解析色散上失误，不再用它做判据。）

**判读矩阵命中行**：「离散色散收敛 + S5b-CW 收敛 -> 脉冲 S5b 的 1.9% 是窗口/色散
测量问题；**E6 的 1% 也应按此重测**」。

**当前结论**：

- 体格式的 evanescent 模式在 CW 稳态下是**二阶收敛**的（`-1.07% -> -0.225%`）；
- 之前所有「evanescent 不收敛」的证据都来自**脉冲 + 单窗口 DFT 的测量偏置**；
- E6 的 1280^2 `+1.150%` 同样需要 **CW 重测**才能定性（深探针区 `z >= 0.2` 拟合，
  避开界面瞬态；或用真空侧两探针测 R，那是 E3 的方法）。


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


### 5.0t **深 node 诊断：内部衰减正确（-0.91%），「快 8 倍」是归一化错误**

算例：`tensor_interface` 几何（`eps2 = diag(4,5,6) eps0`）+ 标量 `sigma = 0.01` + 界面 `z=0`，
`E_ext_grid` 载频脉冲（`zc=-4`、`w=1.5`），480^2（`dt = 2.358654337e-11`）。

**深 node 的材质/M1/M2（z=0.5 与 z=1.0 完全相同）**：

```
eps6: 3.541675125e-11  4.427093906e-11  5.312512688e-11   (= diag(4,5,6) eps0)
sig6: 0.01 0.01 0.01 0 0 0
M1  : 0.9933624128  0.9946863962  0.9955700321   = exp(-sigma dt / eps_ii)  ✓
M2  : 0.6637587238  0.5313603788  0.4429967941   = (1 - M1) / sigma        ✓
```

**结论 1：bulk 用的是材料 2 的 `(4,5,6) eps0`，界面 tau 平均没有漏进内部**
（若漏进去会是 `(2.5,3.0,1.714)`）。

**内部衰减剖面（`ln|Ey|` vs z，峰値）**：

| z | peak abs(Ey) | t_peak |
|---|---|---|
| 0.05 | 2.8273e3 | 7.34 ns |
| 0.10 | 2.7087e3 | 7.71 ns |
| 0.20 | 2.4869e3 | 8.47 ns |
| 0.30 | 2.2843e3 | 9.22 ns |
| 0.50 | 1.9277e3 | 10.71 ns |
| 1.00 | 1.2599e3 | 14.46 ns |

```
slope fit (z >= 0.1): -0.8502 1/m  ->  delta_eff = 1.176 m
theory delta = 1.187 m             ->  -0.91%
peak-time increments match c/sqrt(5) = 1.34e8 m/s (0.05 m per 0.37 ns)  ✓
```

**结论 2：内部 evanescent 衰减与传播速度都正确（-0.91%）**。

**结论 3：「比 delta 预测快 8 倍」是归一化错误，撤回。**
当时用 `ln(1e6 / 1.26e3) = 6.7` 当作 1 m 处的衰减，但 `1e6` 是**源幅度**，不是界面处的场：
`E_ext` 的高斯包络在 z=0 只剩 `1e6*exp(-(4/1.5)^2) ~ 8e2`，而且 `E_ext` 是**全域初值**，
所以 z>0 的探针在脉冲到达前就读到初值的尾（t_peak=7.3 ns < 13.3 ns 的传播时间）。
用 `|Ey(0.5)|` 与 `|Ey(1.0)|` 相除才是纯内部衰减。

**仍开放（与内部衰减无关）**：界面**透射幅度**的绝对定标。要测它需要一个在界面处幅度
已知的入射（例如包络足够宽、源贴近界面，或 CW 稳态 + 单频幅度），属于界面 BC 的定量，
不是 bulk 的问题。


### 5.0u **T1：sigma 张量输入（`sigma_tensor_xx..yz`）已落地**

按 `SIGMA_TENSOR_AND_CONDUCTIVE_INTERFACE.md` §4 的 T1 实现：

| 项 | 内容 |
|---|---|
| H | `m_has_sig_tensor`、`m_sig_tensor[6]`（GpuArray）、`m_sigma_tensor_mf`、6 个 parser |
| 解析 | `sigma_tensor_xx/yy/zz/xy/xz/yz`，逐分量 parser（常数或 `(x,y,z)` 函数） |
| 守卫 | key 存在但无解析 -> abort；全零 -> abort；**与标量 `sigma`/`sigma_function` 同时给出 -> abort** |
| 存储 | cell-centred `m_sigma_tensor_mf`（6 分量，含 `ng_props` ghost），`AllocateLevelMFs` 分配、`InitData` 填充 |

**实测（`xx = 0.01*(1+0.1z)`、`yy=0.02`、`zz=0.03`、`xz=0.004`）**：

```
[MACRO] sigma tensor: xx=0.01 yy=0.02 zz=0.03 xy=0 xz=0.004 yz=0
[SIGT]  z=-1.84375   0.00815625 0.02 0.03 0 0.004 0
        解析 0.01*(1+0.1*(-1.84375)) = 0.00815625   ✓（空间变化正确）
冲突检查: `sigma = 0.1` + `sigma_tensor_*` -> "give either `sigma`/`sigma_function` or ..." abort  ✓
P0 回归 (sigma=0, 512^2): R_flux=0.1335 T_flux=0.8513  逐位不变
```

**实现中抓到的一个真 bug**：tensor 填充块原本被 `if (m_has_eps_tensor || m_has_eps2_tensor)`
守卫，只给 `sigma_tensor_*` 时整块被跳过（σ 张量保持 0）。守卫已改为包含 `m_has_sig_tensor`。

**尚未消费**：`m_sigma_tensor_mf` 目前只是存好，node σ 规则（T2）和 `A = eps^-1 sigma`（T3）
仍用标量 σ。所以 T1 不改变任何仿真结果（P0 逐位不变即为证据）。


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
