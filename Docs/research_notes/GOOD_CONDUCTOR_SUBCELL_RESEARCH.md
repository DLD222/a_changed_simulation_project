# δ ≪ Δx 良导体：文献方法与处理建议（建议稿，非决定）

> 问题：良导体的趋肤深度 δ = sqrt(2/(ω μ σ)) 可能远小于网格 Δx；
> 当前代码（体材料 + 指数步 + subpixel τ）在 δ ≪ Δx 时只退化为 PEC 极限，
> 内部衰减/吸收/相位没有分辨。本文查文献，给出可选方法与建议；**不拍板**。
> 日期：2026-11-04

---

## 0. 结论先行

```
δ ≥ 3Δx（已验收区间）：体材料 + subpixel τ（Oskooi/Deinega 型），二阶；不需要 SIBC
δ ≪ Δx（光学金属/薄导体）：
  首选 SIBC（Leontovich 面阻抗 + 时域实现）
  薄层/薄片：共形 thin-conductor / effective surface impedance
  需要内部场：局部 AMR / subgridding
  色散金属：ADE/Drude/临界点
  subpixel/体平均：只对 δ ≳ Δx 有效
```

---

## 1. 方法分类

### 1.1 SIBC（Leontovich 面阻抗）

频域：
\[
\mathbf E_t = Z_s\,\hat{\mathbf n}\times\mathbf H_t,\qquad
Z_s=\sqrt{\frac{i\omega\mu}{\sigma+i\omega\varepsilon}}
\;\xrightarrow{\text{良导体}}\;
(1+i)\sqrt{\frac{\omega\mu}{2\sigma}}=\frac{1+i}{\sigma\delta}.
\]

- 适用：δ ≪ 曲率半径、δ ≪ λ；
- 时域实现：
  - **Maloney & Smith (1992)**：FDTD 的 SIBC 卷积实现（指数/递归卷积）；
  - **PLRC**（piecewise-linear recursive convolution）：色散/宽带的稳定实现；
  - **ADE / 临界点**：把 \(Z_s(\omega)\) 有理拟合，状态方程推进（宽带、色散金属）；
  - **高阶 SIBC**（Bingle 等 2002；Mitzner）：曲率/二阶修正。
- 优点：网格不必分辨 δ；标准做法（Tidy3D 默认、gprMax 阻抗面、Taflove）；
- 代价：面阻抗近似（Leontovich 一阶；曲率/大角度有误差）；时域卷积/拟合的稳定性与无源性；
- 参考：
  - Maloney & Smith, *IEEE TAP* 1992（[Google Scholar](https://patents.google.com/scholar/16529179475247143740)）
  - Bingle 等, 高阶 SIBC（[IEEE 1229914](https://ieeexplore.ieee.org/document/1229914/references)）
  - 色散层 + PEC 背衬的 SIBC（[Unibo](https://acnpsearch.unibo.it/singlejournalindex/47237365)）

### 1.2 共形薄导体（conformal thin conductor）

- 把薄导体层等效为**面阻抗/薄层**，在共形网格上更新；
- 参考："A conformal scheme for simulation of thin good conducting objects in FDTD"（[IEEE 5171793](https://ieeexplore.ieee.org/document/5171793)）；
- 适用：导体厚度 ≪ Δx 的薄片/涂层；与 ECT/共形几何结合；
- 代价：几何/等效参数的构造；二阶性依赖实现。

### 1.3 Subpixel / 体平均（当前路线）

- **Oskooi & Kottke (2009)**：subpixel smoothing（ε⁻¹ 张量；二阶；[PDF](https://math.mit.edu/~stevenj/papers/OskooiKo09.pdf)）；
- **Deinega & Valuev (2007)**：导电介质的 ε⁻¹ 张量（[Opt. Lett. 32, 3429](https://pubmed.ncbi.nlm.nih.gov/18059956/)）；
- **Meep**：subpixel 只对无损 ε（复 ε 虚部忽略）；
- **Tidy3D**：有损金属默认 SIBC，另有 volumetric averaging；
- 适用：δ ≳ Δx（δ 被 3–6 格分辨，正是本项目已验收区间）；
- 代价：δ ≪ Δx 时无效（退化为 PEC）。

### 1.4 ADE / Drude / 临界点（色散金属）

- 把金属的复 ε 写成 Drude/Lorentz/临界点，ADE 推进；
- 参考：宽带 FDTD（[core.ac.uk](https://core.ac.uk/download/pdf/42812883.pdf)）、统一 FDTD/PML 临界点（[IEEE 6522119](https://ieeexplore.ieee.org/document/6522119)）；
- 适用：色散金属/等离子体；与 SIBC 可组合（SIBC 的色散面阻抗）。

### 1.5 局部 AMR / subgridding

- 局部加密到 δ（或薄层厚度）；需要多级 AMR 或 subgridding；
- 优点：不绕过内部解（符合“不降级”）；
- 代价：AMR 的多级支持（本项目被 `lev==0` 挡住）；subgridding 的稳定性/守恒复杂。

### 1.6 薄线/薄槽 subcell

- 细导线/细缝的 subcell（Holland、thin-strut、Umashankar）；
- 适用：一维细结构；不是一般金属面。

---

## 2. 对比

| 方法 | 适用 | 阶/精度 | 稳定性 | 代价 |
|---|---|---|---|---|
| SIBC（Leontovich） | δ ≪ Δx | 一阶（曲率/大角度有偏） | 被动拟合可稳定 | 卷积/拟合 |
| 高阶 SIBC | δ ≪ Δx，曲率 | 二阶修正 | 同上 | 更复杂 |
| 共形薄导体 | 薄层 | 二阶（实现相关） | 好 | 几何/等效参数 |
| Subpixel/体平均 | δ ≳ 3Δx | 二阶 | 好 | δ≪Δx 无效 |
| ADE/Drude/临界点 | 色散金属 | 二阶（拟合） | 需无源拟合 | 参数提取 |
| 局部 AMR/subgridding | 需要内部场 | 网格阶 | 复杂 | 多级/守恒 |
| 薄线 subcell | 细导线 | 一阶/二阶 | 好 | 仅细结构 |

---

## 3. 对本项目的建议（建议，非决定）

```
1. δ ≥ 3Δx（当前验收）：保持体材料 + subpixel τ；不启用 SIBC（避免绕过内部解）
2. δ ≪ Δx（新需求）：
   - 首选 SIBC：Leontovich + 时域 PLRC/ADE；被动拟合保证稳定；
   - 薄层：共形 thin-conductor / effective surface impedance；
   - 需要内部场：局部 AMR（先把多级 AMR 打开）或 subgridding；
   - 色散金属：ADE/Drude/临界点（可与 SIBC 组合）
3. 验收：
   - 法向入射的复 ε Fresnel（精确）作基准；
   - SIBC 的 δ≪Δx 极限应与 PEC 一致，但吸收/相位要与精确解比；
   - 网格加密收敛阶（SIBC 预期一阶；高阶 SIBC/共形预期二阶）
```

实现落点：

```
EvolveE/MacroscopicEvolveE：在界面/面元上施加 E_t = Z_s n×H_t
时域：PLRC/ADE 的状态量（每个面元），与现有的指数步并存
开关：medium.sibc = on/off；默认 off（δ 被分辨的算例不受影响）
测试：法向入射复 ε Fresnel；良导体极限 R→1/T→0；网格收敛
```

---

## 4. 参考

1. Maloney & Smith, FDTD SIBC (1992)：<https://patents.google.com/scholar/16529179475247143740>
2. 高阶 SIBC（Bingle 等）：<https://ieeexplore.ieee.org/document/1229914/references>
3. 色散层 + PEC 的 SIBC：<https://acnpsearch.unibo.it/singlejournalindex/47237365>
4. 共形薄导体：<https://ieeexplore.ieee.org/document/5171793>
5. Oskooi & Kottke, subpixel smoothing (2009)：<https://math.mit.edu/~stevenj/papers/OskooiKo09.pdf>
6. Deinega & Valuev, 导电介质 ε⁻¹ 张量 (2007)：<https://pubmed.ncbi.nlm.nih.gov/18059956/>
7. 宽带 FDTD：<https://core.ac.uk/download/pdf/42812883.pdf>
8. 统一 FDTD/PML 临界点：<https://ieeexplore.ieee.org/document/6522119>
9. Tidy3D 有损金属：<https://tidy3d-help-center.simulation.cloud/flex_rf_gui/TCM/configuration/advanced/subpixel_averaging/LossyMetal.html>
10. Meep subpixel：<https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/>

---

## 5. 一句话

**δ ≪ Δx 的正解是 SIBC（Leontovich + 时域 PLRC/ADE，曲率大时用高阶 SIBC），薄层用共形 thin-conductor，需要内部场时用局部 AMR/subgridding，色散金属用 ADE/Drude/临界点；subpixel/体平均只对 δ ≳ 3Δx 有效（本项目已验收区间）。建议：δ 被分辨时保持现状（不启用 SIBC），δ ≪ Δx 时新增 SIBC 开关，用法向入射复 ε Fresnel 验收，SIBC 预期一阶、高阶 SIBC/共形预期二阶。**
