# 添加流体求解器并与场求解器耦合：研究与方案（建议稿，非决定）

> 目的：在现有 WarpX（PIC + Yee/PSATD 场求解器 + 冷流体 + 1D 矩测试床）上，
> 添加一个 **5 矩流体求解器**，并与**电磁场求解器**正确耦合。
> 本文是调研与方案建议，**不拍板**；所有“建议/可选”由负责人决定。
> 日期：2026-11-04

---

## 0. 结论先行

```
1. 代码里已有两套“半个答案”：
   - 3D 冷流体（Fluids/）：有电荷/电流沉积、Lorentz 推、Muscl 对流，但无压强/能量；
   - 1D 5 矩（MomentsSolver/）：有 5 矩 + Strang 源，但只在一维、源是简化 Ampere/Lorentz。
2. 场求解器只认两个源：current_fp（J）与 rho_fp（ρ）；冷流体已经接上，5 矩还没接。
3. 耦合的正确顺序是“源半步 → 流体输运 → 场推进 → 源半步”（Strang），
   且 J 必须在半格（n+1/2）落到 Yee 的 E 位置、ρ 落到节点，保证电荷守恒。
4. 显式耦合的硬限制是 ωp·dt（等离子体频率）；超了就要 IMEX/CN 或 hybrid（Ohm 定律）。
```

---

## 1. 现状（代码事实）

### 1.1 场求解器

```
EvolveB(dt, SubcyclingHalf, ...)      // B^{n} → B^{n+1/2} / B^{n+1}
EvolveE(dt, ...) / MacroscopicEvolveE // E^{n} → E^{n+1}
源：current_fp（J，Yee 位置）、rho_fp（ρ，节点）
PSATD：J/ρ 经 Vay 构造、半格对齐
```

### 1.2 3D 冷流体（`Fluids/WarpXFluidContainer`）

```
Evolve():
  DepositCharge(rho_fp, icomp=0)        // 电荷
  GatherAndPush(E_aux,B_aux)            // Lorentz（Higuera–Cary）
  AdvectivePush_Muscl()                 // N, NU 对流
  DepositCharge(rho_fp, icomp=1)
  DepositCurrent(current_fp)            // 电流
```

- 状态：`N`（密度）、`NU`（动量密度）；**无压强/能量**；
- 已在 3D 接上 `rho_fp`/`current_fp`，但 `rho_fp` 只在静电/PSATD 等路径分配（见 §1.4）。

### 1.3 1D 5 矩（`MomentsSolver`）

```
Euler5: U=[rho, rho u, E], p=(γ-1)(E-½ρu²)
MomentCoupling::strang_step:
  source(0.5dt): E_new = E - dt J/eps0 ; u_new = u + dt q/m (E_new + u×B)
  transport(dt): Euler5 (流体) + Maxwell6 (场)
  source(0.5dt)
```

- 有 5 矩、有源、有 Strang；**但只在一维**，场是 6 分量 FV，不是 Yee；
- J 由 `q/m * ρu` 直接算，没有 Yee 位置的沉积/守恒。

### 1.4 缺口

```
- 5 矩没有 3D 实现，也没接 current_fp/rho_fp；
- 冷流体的 rho_fp 在“显式电磁 + 宏观介质”路径未分配（Fluid 的 DepositCharge 被跳过）；
- 电流/电荷的 Yee 沉积（位置、半格、守恒）未做；
- 显式源的 ωp·dt 限制未处理；
- EOS/摩擦/闭合未接（见相关文档）。
```

---

## 2. 物理模型

5 矩（每组分）：
\[
\partial_t\rho+\nabla\cdot(\rho\mathbf u)=0,
\]
\[
\partial_t(\rho\mathbf u)+\nabla\cdot(\rho\mathbf u\mathbf u+p\mathbb I)
=nq(\mathbf E+\mathbf u\times\mathbf B)+\sum_t\mathbf R_{st},
\]
\[
\partial_tE+\nabla\cdot[\mathbf u(E+p)]
=nq\,\mathbf u\cdot\mathbf E+\text{friction work/heat},
\]
\[
p=p(\rho,e),\qquad \mathbf J=qn\mathbf u=\frac{q}{m}\rho\mathbf u,\qquad \rho_c=qn.
\]

场：Maxwell（Yee/PSATD），源为 \(\mathbf J,\rho_c\)。

---

## 3. 耦合方式（核心）

### 3.1 顺序：Strang（建议）

```
source(dt/2):  用 E^n/B^n 推流体动量/能量半步；同时把 J^{n+1/2} 沉积到 current_fp
transport:     流体输运 dt（FV）；场推进 dt（EvolveB/EvolveE）
source(dt/2):  再推半步
```

- 与现有的 `MomentCoupling::strang_step` 一致；
- 场推进内部仍是 `EvolveB(0.5dt, FirstHalf) → EvolveE(dt) → EvolveB(0.5dt, SecondHalf)`；
- J 在 **n+1/2** 落到 E 位置，与 Yee 的 `E^{n+1}=E^n+dt(curlB^{n+1/2}-J^{n+1/2})` 对齐。

### 3.2 沉积：电荷守恒（关键）

```
ρ_c → 节点（rho_fp）
J   → E 位置（current_fp：Jx 在 x 面、Jy 在 y 面、Jz 在 z 面）
```

- 用“流体版本的 Esirkepov/电荷守恒沉积”：由同一组 \((\rho,\mathbf u)\) 构造 J，使
  \[
  \partial_t\rho_c+\nabla\cdot\mathbf J=0
  \]
  在离散意义下成立；
- 否则 Gauss 律漂移，`rho_fp`/`current_fp` 不自洽。

### 3.3 源：Ampere + Lorentz

- 显式（现有 1D 版）：
  \[
  \mathbf E_{\rm new}=\mathbf E-dt\,\mathbf J/\varepsilon_0,\qquad
  \mathbf u_{\rm new}=\mathbf u+\frac{q}{m}dt(\mathbf E_{\rm new}+\mathbf u\times\mathbf B);
  \]
- 场的完整推进由 `EvolveB/EvolveE` 做（含 curl），源只提供 J；
- 能量：源步更新总能量（动能+内能），保持内能非负（现有 1D 版做法）。

### 3.4 稳定性的两条路

| 条件 | 做法 |
|---|---|
| \(\omega_p dt \lesssim 1\) | **显式 Strang**（简单；与 1D 版一致） |
| \(\omega_p dt \gg 1\) | **IMEX/CN**：E 的源项隐式（\(E^{n+1}=E^n-dt J^{n+1/2}/\varepsilon_0\) 的线性化），或 **hybrid**（电子流体 Ohm 定律 + 离子 PIC） |

- 还有 CFL：\(\Delta t\le \min(\Delta x/(c+\dots))\)（场）与 \(\Delta x/(|u|+c_s)\)（流体）；
- 正性：\(\rho>0,\ p>0\)；EOS 的 \(c_s^2>0\)。

---

## 4. 三条实现路线

### 路线 A：3D 5 矩作为新流体/物种（推荐起步）

```
新模块：5 矩状态（ρ,ρu,E）+ FV 3D（MUSCL/HLL）+ EOS + 摩擦
耦合：source(dt/2) 沉积 J/ρ → EvolveB/EvolveE → source(dt/2)
场：现有 Yee（current_fp/rho_fp）
```

- 优点：与现有场求解器/回归兼容；5 矩可独立测试；
- 代价：3D 的 FV、沉积、边界、AMR 要写；
- 适用：从“1D 验证 → 3D 生产”。

### 路线 B：扩展现有 `WarpXFluidContainer` 到 5 矩

```
在冷流体上加 E（能量）/p（压强）/EOS/摩擦；保留现有 GatherAndPush/Deposit
```

- 优点：复用现有 3D 沉积/Lorentz/对流骨架；
- 代价：改旧文件；冷/热两种模型要共存；
- 适用：想最快在 3D 里有“热流体”。

### 路线 C：MomentsSolver 1D → 3D（新驱动/新 app）

```
把 1D FV/Strang 推广到 3D；用公共步骤接现有场求解器
```

- 优点：测试床干净；与 1D 回归一致；
- 代价：3D 的场/边界/诊断要接；
- 适用：以新物理为主、旧特性可后接。

---

## 5. 验证（建议顺序）

```
1. 单组分 Langmuir 振荡：ω = ωp（场+流体耦合的最基本检验）
2. EM 波在等离子体中的色散：ω² = ωp² + c²k²
3. 双流体（电子+离子）：离子声波、双流不稳定性
4. Sod/激波：纯流体 FV 正确性（无场）
5. 电荷/能量守恒：∇·D − ρ、总能量、总动量
6. 与 1D MomentsSolver 的同一算例对比（数值门）
```

---

## 6. 文献与参考

- WarpX 冷流体模型（现有实现）：<https://github.com/blast-warpx/warpx/blob/e0af5f84c3fee8dd332b0b6b274ea22d175fda73/Docs/source/theory/cold_fluid_model.rst>
- Gkeyll 5/10 矩两流体（方程与闭合）：<https://ar5iv.labs.arxiv.org/html/1903.09618>
- 两流体（含 Prandtl 修正/摩擦）：<http://export.arxiv.org/pdf/2005.06295>
- 电子流体 + 压强张量 + 摩擦 \(\mathbf R_{ei}\)：<https://export.arxiv.org/pdf/2503.03387v1>
- 3D 磁约束两流体数值方法：<https://julib.fz-juelich.de/vufind/EdsRecord/edsbas,edsbas.EB1E7E63>
- Hybrid fluid-PIC（Ohm 定律 + 离子 PIC）：<https://scholar.google.com.vn/citations?view_op=view_citation&hl=vi&user=WR67sz8AAAAJ&sortby=pubdate&citation_for_view=WR67sz8AAAAJ:CMvovTBb2okC>
- 电流分解 \(I=qn u+J\)（沉积/守恒参考）：<https://pure.mpg.de/rest/items/item_147994_2/component/file_147993/content?download=true>

---

## 7. 风险与决策点（留给负责人）

| # | 决策点 | 选项 | 备注 |
|---|---|---|---|
| 1 | 5 矩实现位置 | 新模块 / 扩展冷流体 / 1D→3D | A/B/C |
| 2 | 耦合顺序 | Strang / 同时 / IMEX | 显式先做 |
| 3 | 电流沉积 | 直接 `q/m ρu` / 电荷守恒沉积 | 后者保 Gauss |
| 4 | 场求解器 | 现有 Yee / PSATD | PSATD 的 J/ρ 半格对齐 |
| 5 | 刚性 | 显式 / IMEX / hybrid | 由 ωp·dt 定 |
| 6 | AMR/边界 | 先 uniform / 后 AMR | 影响工作量 |
| 7 | EOS/摩擦 | 理想气体 / 外部 EOS / 三模型摩擦 | 见相关文档 |

---

## 8. 一句话

**现有代码已有“3D 冷流体（沉积/Lorentz/对流）”和“1D 5 矩（Strang 源）”两半；把 5 矩接到场求解器的关键是：状态用 \((\rho,\rho\mathbf u,E)\)、源用 \(\mathbf J=qn\mathbf u\) 与 \(\rho_c=qn\)、按 Strang「源半步 → 流体输运 → 场推进 → 源半步」耦合、J 沉积到 E 位置（n+1/2）且电荷守恒、显式受 \(\omega_p dt\) 限制（超了上 IMEX/hybrid）；建议先 1D 验证（Langmuir/EM 色散/双流体），再按 A（新 3D 模块）或 B（扩展现有冷流体）落地。**
