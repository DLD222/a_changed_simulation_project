# σ 张量与导电界面平滑：文献调研、更正与实施路线

> 目的：为"ε/μ/σ 3×3 张量"中仍缺的 **σ 张量输入**与**导电界面处理**给出有文献依据的路线，
> 并更正此前两处不准确的技术判断。所有结论都对应到可核查的来源或源码。

---

## 0. 现状与问题

| 项 | 现状 |
|---|---|
| ε 3×3 张量（parser + cell MF + Kottke τ） | ✅ 已交付 |
| σ 矩阵指数 `M1/M2`（数学层支持 3×3） | ✅ 已交付 |
| **σ 张量输入**（`sigma_tensor_*`） | ❌ 只有标量 `m_sigma` |
| **μ 3×3 张量**（`H=B/μ` 的 μ⁻¹ 行） | ❌ 只有标量 |
| σ 的界面处理 | 当前是 **staircase**（切向 0 / 法向调和） |

问题：σ 的界面该用 staircase，还是像 ε 一样的 τ/ε⁻¹ 张量平滑？此前曾判断"real-σ 的 τ 没有二阶依据、会带回自锁"，本文核查后**更正**。

---

## 1. 文献结论

### 1.1 Deinega & Valuev (2007)：ε⁻¹ 张量法**明确适用于导电介质**

> A. Deinega, I. Valuev, "Subpixel smoothing for conductive and dispersive media in the
> finite-difference time-domain method," *Opt. Lett.* **32**, 3429 (2007).

摘要原文：

> "There are different approaches for solving this problem, and the most successful of them
> are based on correct approximation of **inverse permittivity tensor ε⁻¹** at the material
> interface. We report an application of this tensor method for **conductive and dispersive
> media**."

来源：[PubMed 18059956](https://pubmed.ncbi.nlm.nih.gov/18059956/)、[Optica OL 32, 3429](https://opg.optica.org/abstract.cfm?URI=ol-32-23-3429)。

**含义**：对导电界面做"复 ε 的 ε⁻¹ 张量（即 τ）平滑"是**有文献依据的二阶路线**；staircase 是一阶近似。此前"real-σ 的 τ 没有依据"的说法**不成立**。

### 1.2 Meep：子网格平滑只对**无损** ε

> "It only applies to frequency-independent, **lossless** dielectrics... For materials with a
> complex-refractive index, the imaginary part is simply ignored."

来源：[Meep Subpixel Smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/)。

**含义**：Meep 的 τ 公式（`ε̃⁻¹ = P⟨ε⁻¹⟩ + (I−P)⟨ε⟩⁻¹`）本身不覆盖 σ；但它给出的是 **ε⁻¹ 张量**的形式，Deinega 把同一形式推广到复 ε。

### 1.3 Tidy3D：有损金属的三种处理

| 选项 | 说明 |
|---|---|
| **SurfaceImpedance（默认）** | 对 `LossyMetalMedium` 施加一阶表面阻抗 BC |
| Staircasing | 阶梯化 |
| **VolumetricAveraging** | 对 Yee 单元体积平均材料属性 |

来源：[Tidy3D Lossy Metal Subpixel Averaging](https://tidy3d-help-center.simulation.cloud/flex_rf_gui/TCM/configuration/advanced/subpixel_averaging/LossyMetal.html)。

**含义**：工业界对"δ 分辨不了"的有损金属默认 SIBC；**体积平均**是可选路线。我们的 δ=9.2 cm、dx=3.1/1.6 cm（3–6 格/δ）是**被分辨**的，所以既不需要 SIBC（会绕过内部解，降级），也应该用**张量平滑**而不是 staircase。

### 1.4 Oskooi/Kottke (2009)：ε⁻¹ 张量的二阶平滑理论

> A. Oskooi, A. Kottke, S. Johnson, "Accurate finite-difference time-domain simulation of
> anisotropic material by subpixel smoothing," *Opt. Lett.* **34**, 2778 (2009).

来源：[MIT PDF](https://math.mit.edu/~stevenj/papers/OskooiKo09.pdf)。给出无损各向异性 ε 的 ε⁻¹ 张量平滑及其二阶精度；是 Deinega 复 ε 推广的基础。

### 1.5 教科书与共形方法

- A. Taflove, S. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain
  Method*：材料界面的平均（法向调和、切向算术）、staircase 误差、表面阻抗 BC、共形技术。
- Dey & Mittra (1997)：局部共形 FDTD（PEC 边界的一阶修正）；属于"边界几何修正"，与复 ε 平滑互补。

### 1.6 小结

| 路线 | 依据 | 适用 |
|---|---|---|
| **复 ε 的 ε⁻¹ 张量（τ）平滑** | Deinega 2007；Oskooi 2009（无损理论） | **δ 被分辨**（我们的情况），二阶 |
| staircase | 当前实现 | 一阶，最简单 |
| SIBC（Leontovich） | Tidy3D 默认；Taflove | δ ≪ dx（光学金属） |
| 共形 FDTD | Dey–Mittra | PEC 边界几何 |

---

## 2. 两处更正

### 2.1 "real-σ 的 τ 没有二阶依据" → **错**

Deinega 2007 明确把 ε⁻¹ 张量法用于导电/色散介质。正确表述是：**对复 ε 做 ε⁻¹ 张量平滑**（等价于对 ε 与 σ 一起平滑），而不是"对实 σ 单独做 τ"。

### 2.2 "指数步的 J 必须是 S·E" → **错**

代码中指数步的 `rhs = curlH − J` 里，**J 是外部/粒子电流**（`current_fp`）；**传导电流 σE 已经在 `A = ε⁻¹σ` 的 `M1/M2` 里**。

```
正确：A = ε⁻¹σ（3×3），M1 = exp(−A dt)，M2 = (1−e^{−A dt})A⁻¹ε⁻¹；J 保持外部电流
错误：把 S·E 再加进 rhs（传导电流会被算两次）
```

`S·E` 只用于 **ρ_medium 的电荷更新**：`ρ += −dt ∇·(σE)`。

### 2.3 自锁归因更正

自锁根因是 **BE 在界面 node 用 cell 插值 σ（≈σ/2）过阻尼 + 节点副本/grid 边界 bug**（已由"BE 用 node σ"和"逐分量 node 指数步"修复，V1 逐位）。**τ 阻尼本身不是自锁原因**——在修复后，界面 node 的 σ_eff>0 走指数步，其 E_t 衰减是**物理表面损耗**。

---

## 3. 正确的 formulation

### 3.1 复 ε 的界面平滑（Deinega 路线）

对每个 node（或 cell，取决于用于哪条路径）：

```
ε̃ = ε + i σ / ω          （符号需与代码 e^{−iωt} 约定核对；亦可写成 ε − iσ/ω）
ε̃_node = τ^{-1}( ⟨τ(ε̃)⟩ )  （沿界面法向 n 的 ε⁻¹ 张量平滑）
ε_node = Re(ε̃_node)
σ_node = ω · Im(ε̃_node)    （张量，6 分量）
```

- **bulk（无界面）**：ε_node = 局部 ε 张量，σ_node = 局部 σ 张量；
- **界面**：用上面的复 ε τ；这同时平滑了 ε 和 σ（切向/法向自动按 P 投影处理），**不需要**单独的 staircase。

### 3.2 指数步

```
A = ε_node^{-1} σ_node            （3×3；ε_node SPD）
M1 = exp(−A dt) = ε^{-1/2} exp(−B dt) ε^{1/2},  B = ε^{-1/2} σ ε^{-1/2}
M2 = ε^{-1/2} (1−e^{−B dt}) B^{-1} ε^{-1/2}
E^{n+1} = M1 E^n + M2 (curlH − J)      J = 外部/粒子电流
```

### 3.3 电荷步

```
J_cond = σ_node · E                （3×3 × 3-vector）
ρ_medium += −dt ∇_h · J_cond
```

这是"σ 张量 + Gauss 记账"的后续步骤；先做 bulk 验证，再接电荷。

---

## 4. 实施路线

| 步 | 内容 | 验收 |
|---|---|---|
| **T1** | `macroscopic.sigma_tensor_xx..yz`（6 parser）+ cell MF（6 分量，含 ghost）；与标量 `sigma` 冲突检查、全零 abort | 解析值正确 |
| **T2** | 复 ε 的 node τ（界面）/局部张量（bulk）→ `ε_node`、`σ_node`（6 分量） | 标量退化：`diag(算术,算术,调和)`；E3 `diag(2.5,3.0,1.714)` |
| **T3** | `A = ε_node⁻¹σ_node`；`M1/M2` 用现有 `sigma_exponential_matrices`；J 不动 | vs Python 3×3 expm ≤1e-10 |
| **T4** | bulk 各向异性 σ 验证：对角 → 各主轴 `exp(−σ_i t/ε_i)`；旋转 → 本征模 | 解析/本征值 |
| **T5** | 界面 + σ 张量：R/T 与"复 ε 的 lossy Fresnel"比 | ≤3%（32 ppw） |
| **T6** | `ρ_medium` 的张量电流 `σ_node·E`；Gauss 一致性 | 残差 `4e-25`（内部） |
| **T7** | 回归：σ=0 逐位、S5 `−0.215%`、E3、E6 | 逐位/数值 |

---

## 5. 现状与影响

| 项 | 结果 |
|---|---|
| 体格式（无界面 evanescent） | S5b-CW `−1.070% → −0.225%`（resid `5e-4`）✅ |
| 界面反射 | `R_flux = 0.3197` vs 精确 `0.3059`（+4.5%）✅ |
| E6 斜率 | 深区 640² `−1.58%`、1280² `−0.61%`（量级符合二阶；拟合区间敏感 ±1%） |
| 自锁 | 已修复（BE node σ + 逐分量 node 指数步，V1 逐位） |
| **σ 张量 / μ 张量** | **未交付**（本文给出 σ 的路线；μ 需 μ⁻¹ 行） |

---

## 6. 风险与约定

1. **符号约定**：`ε̃ = ε ± iσ/ω` 取决于代码的 `e^{±iωt}`；实现前先用一个均匀算例核对 `Im(ε̃)` 的符号与衰减方向。
2. **界面帧**：复 ε 的 τ 需要界面法向 `n = ∇φ/|∇φ|`；与 ε 的 Kottke 使用同一 `frame/tau/tau_inv`。
3. **bulk vs 界面**：bulk 用局部张量（无 τ）；只有界面 node 用复 ε τ。
4. **Gauss**：张量 σ 的电荷步必须与 ε_node 一致；否则 `∇·D − ρ` 会漂移。
5. **不降级**：不用 SIBC（δ 被分辨）；不用 staircase 替代二阶 τ（staircase 仅作 fallback）；J 不与传导电流重复计算。

---

## 7. 参考

1. A. Deinega, I. Valuev, *Opt. Lett.* **32**, 3429 (2007): <https://pubmed.ncbi.nlm.nih.gov/18059956/>
2. Meep, *Subpixel Smoothing*: <https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/>
3. Tidy3D, *Lossy Metal Subpixel Averaging*: <https://tidy3d-help-center.simulation.cloud/flex_rf_gui/TCM/configuration/advanced/subpixel_averaging/LossyMetal.html>
4. A. Oskooi, A. Kottke, S. Johnson, *Opt. Lett.* **34**, 2778 (2009): <https://math.mit.edu/~stevenj/papers/OskooiKo09.pdf>
5. A. Taflove, S. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain Method*.
6. K. Dey, R. Mittra, "A locally conformal finite-difference time-domain algorithm for modeling perfectly conducting objects," *IEEE Microw. Guided Wave Lett.* **7**, 273 (1997).
7. A. Kottke, *Phys. Rev. E* **77**, 036611 (2008).
