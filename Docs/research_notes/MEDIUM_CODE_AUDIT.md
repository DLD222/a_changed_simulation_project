# 源码介质（dielectric/material）部分审计：问题、病根与修复路线

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 日期：2026-10-31 · 范围：本项目 fork 的介质实现（`MacroscopicProperties`、`MacroscopicEvolveE.cpp`、
> `EvolveECartesianECT` 的 M2b-1 材料分支、`FieldAccessorMacroscopic`、`Interp`、测试与文档）
> 方法：通读源码 + 与上游 WarpX 对比 + 对照文献（Taflove 2005、Hwang–Cangellaris 2001、Kottke 2008、
> Meep、gprMax）+ 本项目已有 Mie/TE-TM 数据。**本审计只读代码，未修改。**
>
> 结论：**三条病根**——
> **A. 介质系数是"cell-centered + 算术平均"的一阶处理**（法向分量应调和/张量，界面反射是一阶误差）；
> **B. μ 在两个路径中处理不一致**（Yee 路径 `H=B/μ` 有半格位置错误；ECT 材料路径**完全忽略 μ**）；
> **C. 介质路径没有闭环**（EB 用阶梯掩膜、无界面回归测试、PML 与介质不匹配、文档未声明限制）。
> 当前测试只覆盖"均匀 εr=1.5"，所以 A/B 全部测不出来。

---

## 0. 问题清单（按严重性）

| 编号 | 问题 | 源码位置 | 病根 | 严重性 |
|---|---|---|---|---|
| P1 | **ECT 材料路径忽略 μ**（`dE` 只用 `c²∇×B/εr`） | `EvolveE.cpp:377-387`（Ex）、`:400-409`、`:425-434` | B | **高**（μr≠1 静默错） |
| P2 | `H=B/μ` 为**逐点相除**，Bx/Bz 处 μ 取在格心（半格偏移） | `FieldAccessorFunctors.H:48`；调用 `MacroscopicEvolveE.cpp:168-170` | B | **高**（TE 的 Bx/Bz 受影响） |
| P3 | ε/σ 一律**算术平均**到 E 位置；法向分量未用调和/张量 | `Interp`（`ablastr/coarsen/sample.H:49+`）；调用 `MacroscopicEvolveE.cpp:186-190` | A | **高**（界面一阶） |
| P4 | `algo.macroscopic_sigma_method` 在 ECT 路径**被忽略**（恒用 CN） | `MacroscopicEvolveE.cpp:61-64` → `EvolveECartesianECT` | C | 中 |
| P5 | 调试打印 `[DIAG-ECT]`/`[DIAG-S]` + `static` 计数器在并行区 | `EvolveE.cpp:288-297`、`:352-364` | C | 中（竞态+刷屏） |
| P6 | **无介质界面回归测试**（仅均匀 εr=1.5） | `Examples/Tests/macroscopic_solver/`、`embedded_boundary_cube/*_macroscopic` | C | **高**（掩盖 P1–P3） |
| P7 | assert 文案写"relative permittivity"，实际单位是**绝对 ε** | `MacroscopicProperties.cpp:159-161`；示例 `1.5*epsilon0` | C | 中（用户陷阱） |
| P8 | PML 与介质不匹配（PML 按真空推导，未检查/未声明） | `WarpXEvolve.cpp:660-668`（`DampPML` 无条件） | C | 中 |
| P9 | Yee 材料路径 + EB 是**阶梯**（无 ECT 双对偶几何） | `MacroscopicEvolveE.cpp` 的 Yee 分支 | C | 高（M2b 需求） |
| P10 | ECT 材料分支的 ε 守卫点 `eps_arr(i,j,k)` 与插值点不一致 | `EvolveE.cpp:381`、`:404`、`:429` | A（次要） | 低 |

---

## 1. 介质代码路径总览

```
algo.em_solver_medium = macroscopic
  └─ WarpXEvolve.cpp:649  → MacroscopicEvolveE(dt)
       ├─ m_fdtd_algo == ECT && ect_conformal_E > 0 → EvolveECartesianECT(..., macroscopic_properties)
       └─ 否则（Yee 等）→ MacroscopicEvolveECartesian<T_Algo, T_MacroAlgo>：
            LaxWendroffAlgo / BackwardEulerAlgo（algo.macroscopic_sigma_method 选择）
```

- 材料场：`MacroscopicProperties`：σ/ε/μ **各一个 cell-centered MultiFab**（`MacroscopicProperties.cpp:126-131`），
  常数或 parser 初始化（`:33-115`），`lev==0` 限定（`WarpX.cpp:2660-2662`，无 AMR）；
- Yee 路径：`H = B/μ`（`FieldAccessorMacroscopic`）+ 算术 `Interp` 到 E 位置 + α/β；
- ECT 材料路径（fork 新增，M2b-1）：`EvolveECartesianECT` 内 `inv_epsr = ε0/ε`、CN 的 α/(1+s)，
  **无 μ**；带双对偶几何（`A_y/Lx/Lz`）与 `A_stab` 下限。

---

## 2. 逐项分析

### P1（高）ECT 材料路径忽略 μ

- 源码：`EvolveE.cpp:377` `dE = -c2*dt*(By(i,j)-By(i,j-1))/Ld - c2*dt*mu0*jx`；
  材料更新 `:386-387` `inv_epsr = epsilon0/eps_abs; Ex = alpha_*Ex + inv_epsr*dE/(1+s_)`。
- 推导：正确形式 `∂E/∂t = (1/ε)(∇×B/μ − J)`；代码只乘 `ε0/ε`，**缺 `1/μr`**：
  `c²∇×B/εr` 应为 `c²∇×B/(εr μr)`。
- 后果：`μr≠1` 时静默错误（E 更新偏大 μr 倍），无 assert/无提示；Yee 路径至少用了 `H=B/μ`。
- 文献：Taflove & Hagness 2005（μ 的处理与 ε 对称，`H=B/μ` 必须在 E 更新前算好）；
  Xiao–Liu ECT（共形与材料需同时保持 `D=εE`、`H=B/μ` 的一致性）。
- 修复：ECT 材料路径加 `mu`：`Interp(mu_arr, mu_stag, E_stag)` 或 `H = B/μ` 的边界插值，
  `dE` 用 `c²∇×(B/μr)`；并加 `μr` 扫描回归。

### P2（高）`H=B/μ` 逐点相除（半格偏移）

- 源码：`FieldAccessorFunctors.H:48` `return m_field(i,j,k,ncomp) / m_parameter(i,j,k);`；
  注释却写"mu is cell-centered and will be interpolated/averaged"（`MacroscopicEvolveE.cpp:165-167`）。
- 2D XZ 位置：`Bx=(1,0)`（z 边）、`Bz=(0,1)`（x 边）、`By=(0,0)`（格心）；μ 在格心 `(0,0)`。
  于是 `Bx(i,j)/μ(i,j)` 的 μ 相对 Bx 偏了 `(dx/2, 0)`，`Bz` 偏 `(0, dz/2)` → **一阶位置误差**；
  `By` 恰好正确。
- 后果：μ 变化时 Bx/Bz（TE 的两个分量）的 H 错位；均匀 μ 时无影响（所以现有测试测不出）。
- 文献：Taflove 2005（H 在边上的 μ 平均与 ε 在边上的平均对偶）；Hwang–Cangellaris 2001（法向/切向规则）。
- 修复：`FieldAccessorMacroscopic` 按 B 的 IntVect 取 μ 的 2/4 点平均（与 `Interp` 同规则），或
  预先把 μ 插值到 B 位置存一份。

### P3（高）ε/σ 的算术平均（界面一阶）

- 源码：`Interp` 对 cell-centered → E 交错位置做**算术平均**（`ablastr/coarsen/sample.H:49+`；
  调用 `MacroscopicEvolveE.cpp:186-190/206-210/227-231`，ECT 路径 `EvolveE.cpp:381-383/404-406/429-431`）。
- 理论：在介质界面处，**法向 E 分量**的等效 ε 应取**调和平均**（或 Kottke 的各向异性张量），
  切向才用算术；纯算术平均的界面反射/透射是一阶误差（Hwang–Cangellaris 2001；
  Kottke 2008；Taflove 2005 的界面平均讨论）。
- 后果：Fresnel R/T 在界面处偏离解析值 O(h)，且对大 ε 对比（如 1:16）最严重；现有测试均匀介质测不出。
- 修复（对齐 `M2B_IMPLEMENTATION_RESEARCH.md` §4.2）：
  - 近场：ε̃ 张量（Kottke）或法向调和/切向算术；
  - 结合 EB：用与 `Lx/Lz/A_y` **同源**的 level-set 构造 ε̃；
  - 加 Fresnel 回归（M2b V1–V8）。

### P4（中）`macroscopic_sigma_method` 在 ECT 路径无效

- 源码：`MacroscopicEvolveE.cpp:61-64` 直接调 `EvolveECartesianECT`，不传 `macroscopic_solver_algo`；
  ECT 内核固定 CN（`(1−s)/(1+s)`，`EvolveE.cpp:384-386`）。
- 后果：用户选 `backwardeuler`/`laxwendroff` 在 ECT 路径无效果；文档未说明。
- 修复：ECT 内核接受 `MacroscopicSolverAlgo`，或文档明确"ECT 路径固定 CN"。

### P5（中）调试残留

- `EvolveE.cpp:288-297`：`static int dbg_n` + `amrex::Print()`（在 MFIter 循环内，OpenMP 下竞态）；
  `:352-364`：`static int dbg_s` + `ReduceOps`（在并行区内）→ 重复/交错打印，污染日志。
- 修复：删除（或 `#ifdef` + 单线程守卫）。

### P6（高）无界面回归测试

- 现有：`Examples/Tests/macroscopic_solver/inputs_test_{1d,2d}_macroscopic_solver`（**均匀** εr=1.5 的模式频率）、
  `embedded_boundary_cube/inputs_test_3d_embedded_boundary_cube_macroscopic`（**均匀** εr=1.5）。
- 后果：P1（μ）、P2（μ 偏移）、P3（界面平均）都测不出；M2b 的 Fresnel 验收完全缺失。
- 修复：加"介质平板 R/T vs Fresnel"回归（垂直/斜 s-p、ε 扫描、μ 扫描），阈值 ≤1%。

### P7（中）单位/文案陷阱

- `MacroscopicProperties.cpp:159-161` 的 assert 文案："relative permittivity `epsilon`"，
  但实际输入是**绝对 ε**（示例 `macroscopic.epsilon = 1.5*epsilon0`）。
- 用户若填 2.25（相对值）→ `inv_epsr = ε0/2.25 ≈ 3.9e-12` → E 几乎冻结（静默错）。
- 修复：改文案为 "absolute permittivity"，并在 `ReaddParameters` 里对 `ε < 1e-6*ε0` 之类给警告/断言。

### P8（中）PML 与介质不匹配

- `WarpXEvolve.cpp:660-668`：`DampPML()` 无条件调用；PML 按真空波阻抗推导，
  介质中 `Z=√(μ/ε)≠Z0` → 反射（Oskooi 2008；Taflove 2005 Ch.7）。
- 修复：文档声明"PML 外必须真空"或实现材料一致的 PML（C-PML with material）；至少加断言/警告。

### P9（高，设计缺口）Yee 材料路径 + EB 是阶梯

- Yee 材料分支只有 `eb_update_E` 掩膜（阶梯）；M2b 要求"曲面介质 + 共形"，
  需要把材料并入 ECT 双对偶几何（`EvolveECartesianECT` 已有材料分支，但见 P1–P5）。
- 修复：以 ECT 材料分支为唯一介质路径（修复 P1–P5 后），Yee 材料路径仅保留无 EB 场景。

### P10（低）ε 守卫点不一致

- `EvolveE.cpp:381` `if (has_macro && eps_arr(i,j,k) > 0.0)` 用格心守卫，而更新用插值 ε；
  两者在界面处可能不同。修复：守卫用插值值（已在 `eps_abs` 里 `max(...,1e-30)`）。

---

## 3. 病根归纳

1. **病根 A：介面系数的一阶处理**——cell-centered + 算术平均，未按"法向调和/切向算术"或 Kottke 张量；
   界面反射一阶，Fresnel 验证会直接暴露。
2. **病根 B：μ 的两套处理**——Yee 路径的 `H=B/μ` 半格错位，ECT 材料路径完全忽略 μ；
   `μr` 扫描一测就露。
3. **病根 C：没有闭环**——EB 阶梯、无界面测试、PML 不匹配、文档/单位不清、调试残留；
   即使 A/B 修好，也会被 C 掩盖或放大。

---

## 4. 修复路线（与 `M2B_IMPLEMENTATION_RESEARCH.md` 对齐）

| 步 | 内容 | 涉及 | 验收 |
|---|---|---|---|
| R1 | 删调试残留；ECT 材料路径加 μ；修 `H=B/μ` 位置 | P1/P2/P5 | μr 扫描（1,2,4）R/T 与解析 ≤1% |
| R2 | ε/σ 平均改为法向调和 + Kottke 张量（先 2D isotropic） | P3 | Fresnel 垂直/斜 s-p ≤1% |
| R3 | 材料与 EB 共形统一（ECT 路径 + 同源 level set ε̃） | P9/P3 | 曲面介质 vs Meep ≤2% |
| R4 | 加回归：介质平板 R/T、μ 扫描、界面能量、PML+介质 | P6/P8 | 全套 CI |
| R5 | 文档：单位（绝对 ε）、ECT 路径固定 CN、PML 限制 | P4/P7/P8 | 文档评审 |

---

## 5. 验收数字（目标）

- Fresnel R/T：垂直/30°/45°/Brewster 角，\|R\|、\|T\| 误差 ≤1%，相位 ≤2°；
- ε 扫描 1.5/2.25/4/16；μ 扫描 1/2/4；σ 扫描 0/0.1/1（CN 稳定）；
- 界面能量：无损介质总能量漂移 ≤1e-6（216 步）；
- 曲面介质（圆柱/透镜）与 Meep 对照 ≤2%（同网格/同 level set）。

---

## 6. 参考文献

1. A. Taflove, S. Hagness, *Computational Electrodynamics: The FDTD Method*, 3rd ed., 2005
   （Yee 算法/材料界面平均/μ 的处理；PML Ch.7；共形 Ch.8）。
2. K.-P. Hwang, A. Cangellaris, "Effective permittivities for second-order accurate FDTD equations at
   dielectric interfaces," IEEE MWCL 11(4):158–160, 2001（法向调和/切向算术，二阶）。
3. J. Kottke, A. Farjadpour, S. Johnson, "Perturbation theory for anisotropic dielectric interfaces…,"
   Phys. Rev. E 77, 036611, 2008（各向异性 ε̃ 张量与 σ 处理）。
4. A. Farjadpour et al., Opt. Lett. 31:2972, 2006（subpixel）。
5. A. Oskooi et al., J. Comput. Phys. 227:5971, 2008（PML 的失效与材料）。
6. Meep：`doc/docs/Subpixel_Smoothing.md`（界面场一阶的限制）、`src/anisotropic_averaging.cpp`。
7. gprMax（本机）：阶梯 + 色散材料实现（对照）。
8. 本项目：`M2B_IMPLEMENTATION_RESEARCH.md`（M2b 方法）、`MIE_FIX_METHODS_RESEARCH.md`（测量/几何）。

## 附录：关键行号

| 主题 | 位置 |
|---|---|
| σ/ε/μ cell-centered 分配 | `MacroscopicProperties.cpp:126-131` |
| 常数/parser 初始化 | `MacroscopicProperties.cpp:33-115`、`:137-180` |
| ε>0 assert（"relative" 文案） | `MacroscopicProperties.cpp:159-161` |
| Yee 材料更新（α/β、Interp） | `MacroscopicEvolveE.cpp:180-232` |
| H=B/μ 逐点 | `FieldAccessorFunctors.H:48`；`MacroscopicEvolveE.cpp:168-170` |
| ECT 材料分支 | `MacroscopicEvolveE.cpp:61-64` → `EvolveE.cpp:244-440` |
| ECT 材料更新（ε 因子） | `EvolveE.cpp:377-387`、`:400-409`、`:425-434` |
| 调试残留 | `EvolveE.cpp:288-297`、`:352-364` |
| 无 AMR | `WarpX.cpp:2660-2662` |
| ECT 需 EB | `WarpX.cpp:588` |
| 均匀介质测试 | `Examples/Tests/macroscopic_solver/`、`embedded_boundary_cube/*_macroscopic` |
| PML 无条件 | `WarpXEvolve.cpp:660-668` |

*本审计未修改代码；修复需按 §4 的 R1–R5 逐步提交并加回归。*
</content>
