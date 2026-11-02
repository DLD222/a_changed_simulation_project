# 电导率/介电常数模块：分级算例与物理正确性验证方案

> 目标：对已交付的 ε/σ 模块（标量/张量、界面、subpixel、面电荷）做**由易到难、定量**的验证；
> 每个算例都给出解析参考、运行方法、分析脚本、判据和预期数值。
> 验证顺序按"信息量/成本"排列：先数学核，再体格式，再单界面，再张量，再界面+张量，再电荷，最后斜入射与鲁棒性。

---

## 0. 能力矩阵（验证范围）

| 能力 | 状态 | 对应算例 |
|---|---|---|
| ε 标量/各向同性 | ✅ | L1、L2 |
| ε 张量（对角/旋转） | ✅ | L3、L5 |
| σ 标量 | ✅ | L1、L2、L6 |
| σ 张量（对角/旋转） | ✅（S1/S2） | L4、L5 |
| 界面 subpixel（cell Kottke / node τ） | ✅（默认 node τ，A/B 开关保留） | L3 |
| 有损界面 Fresnel R/T | ✅ | L2、L5 |
| 面电荷 ρ_medium / q_s | ✅（诊断量） | L6 |
| 指数积分 M1/M2、C、J0/J1 | ✅ | L0 |
| 斜入射（p-pol，E_n≠0） | ⚠️ 待验证 | L7 |
| 粒子 + 介质（ρ_p+ρ_medium） | ⚠️ 未接线 | L9（限制） |
| PML + 介质、>2 材料、AMR、μ 张量、σ 复 ε τ | ❌ 未交付 | L10（不测） |

---

## 1. 验证方法论（文献与社区做法）

1. **解析平面波/界面基准**：FDTD 验证的标准做法是"均匀有损介质的衰减+相位"与"单界面 Fresnel R/T"，并与解析式比；见 Taflove & Hagness 的验证章节、NASA TM-2001-210862 的有损介质算例。
2. **收敛阶**：同一算例在 16/32 ppw（或 480²/864²）下误差应按二阶下降；只报单网格的结果不算验证。
3. **能量账**：有损界面用 `R + T + A = 1`（A 由板内吸收解析给出）；不能只用 `R+T=1`。
4. **subpixel 的社区现状**：Meep 的 τ 只对无损 ε（复 ε 虚部忽略，[Meep Subpixel Smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/)）；Tidy3D 对有损金属默认 SIBC、另有 VolumetricAveraging（[Tidy3D](https://tidy3d-help-center.simulation.cloud/flex_rf_gui/TCM/configuration/advanced/subpixel_averaging/LossyMetal.html)）；导电介质的 ε⁻¹ 张量方法见 Deinega & Valuev（[Opt. Lett. 32, 3429 (2007)](https://pubmed.ncbi.nlm.nih.gov/18059956/)）。
5. **不降级原则**：δ 被分辨（δ/dx ≥ 3）时不用 SIBC；staircase 只作 fallback；复 ε τ 是可选升级。
6. **测量陷阱（本项目已踩过，必须写进流程）**：
   - `macroscopic.sigma` 是**全域单场**，要"只在某材料"必须用 `sigma_function`；
   - 判断界面反射前，先做**路径预算**（源→界面的 e^{−∫dz/δ}），入射到不了界面 = 算例作废；
   - 真空参考跑**没有 slab 吸收**，其 PEC 回波比有损跑早得多，T 的窗口必须停在参考回波之前；
   - 定点 dump 的 index 必须在 FAB 的 **valid box** 内；
   - 矩阵验证的 `dt` 必须取 `run.log` 的**实际值**（手算近似会带来 1e-4 级偏差）；
   - 绝对幅度不能拿源幅度归一（界面处的包络远小于源幅度）。

---

## 2. 分级算例

### L0 — 数学核（无仿真，秒级）

**目的**：把指数积分与电荷平均核钉死，避免"物理对但矩阵错"。

| 子项 | 内容 | 判据 |
|---|---|---|
| L0.1 | `M1/M2` vs Python 3×3 `expm`（ε 标量/张量、σ 标量/张量、旋转） | ≤1e-10 |
| L0.2 | `C = (1−e^{−Adt})(Adt)⁻¹` vs Python | ≤1e-10 |
| L0.3 | `J0=σC`、`J1=ε(I−C)ε⁻¹` vs Python | ≤1e-10 |
| L0.4 | 极限：σ=0 → `M1=I, M2=dt ε⁻¹`；σ→∞ → `M1→0, M2→1/σ` | 逐位/解析 |
| L0.5 | `J1` 与 `I−C` 的差（张量 ε 下应非零，验证系数纠正） | 记录差值（~2.3e-3） |

**运行**：`benchmarks/cases/c03_medium/math/` 下的 Python 脚本（`python3 test_reference_math.py`），C++ 侧用 `WARPX_DEBUG_NODE` 打印节点矩阵对比。

---

### L1 — 均匀有损介质（无界面）

**物理**：`k = ω√(με)·√(1 − iσ/(ωε))`，`δ = 1/Im k`，`η = √(μ/ε)/√(1 − iσ/(ωε))`。

**配置**：

```
2D XZ，均匀 ε、σ、μ；无界面（不设 interface_function）
源：E_ext 单向脉冲（Bx=-Ey/c，载频 ω）或 CW
探针：沿 z 两点，间距 ≥ λ/2（测相位）且 ≤ δ（测衰减）
网格：16 ppw 与 32 ppw 两档
```

**分析**：

```
slope = d ln|Ey|/dz → δ_eff = 1/|slope|，与 δ 比
phase = Δφ/Δz → v_φ = ω/(dφ/dz)，与 1/Re k 比
```

**判据**：`δ_eff` 误差 ≤2%（32 ppw）、≤5%（16 ppw）；相位 ≤1°；收敛近二阶。

**已测基线**：ε=5ε0、σ=0.01 → `δ_eff=1.176` vs `1.187`（−0.91%）；`v=c/√5`。

---

### L2 — 单界面、标量 ε/σ、法向入射（s-pol）

**物理**（半空间）：

```
η_i = √(μ_i/ε_i)/√(1 − iσ_i/(ωε_i))
r = (η2 − η1)/(η2 + η1),  R = |r|²,  T = 1 − R（半空间）
板：T = |t12 t21 e^{ik2 d}/(1 − r21² e^{2ik2 d})|²,  A = 1 − R − T
```

**配置**：

```
ε1=ε0（真空），ε2=4ε0 或 diag(4,5,6)ε0；界面 z=0
sigma_function = "0.01*(z>0)"（σ 只在材料 2！不能用全局 sigma）
源：E_ext 单向脉冲 + 载频；真空参考跑同几何
探针：z=−2（反射）、z=0.05（透射，导体内部）
窗口：反射窗 [14,24] ns；透射窗必须停在参考跑 PEC 回波之前
```

**判据**：`R` ≤3%（理论 0.1479）；`T` ≤3%（理论 0.852）；`R+T+A=1` ±3%。

**已测基线**：`R=0.1552 (+4.9%)`、`T=0.870 (+2.1%)`、`R+T+A=1.049`（共同归一化偏置 ~5%）。

---

### L3 — 张量 ε（对角/旋转）、σ=0

**物理**：s-pol 只感知 `ε_yy`；`R/T` 用 `η_yy` 的 Fresnel；旋转 ε 时 `ε_yy` 会变，需重算解析。

**配置**：

```
ε1=ε0；ε2=diag(4,5,6)ε0（或绕 x/y/z 旋转 30° 的版本）
σ=0；界面 z=0；s-pol
16/32 ppw 两档
```

**判据**：32 ppw `R` ≤1%、`T` ≤1%；两分支 ε 等价测试（σ=0 vs σ=1e-30）≤1e-12。

**已测基线**（32 ppw）：cell Kottke `+0.86%/+0.24%`；node τ `−0.86%/+0.53%`（默认 node τ）。

---

### L4 — 张量 σ（对角/旋转）、bulk

**物理**：无界面时 `σ_node = 局部张量`（精确）；模式的衰减率是 `A=ε⁻¹σ` 的本征值 `λ_i`（或 `B=ε^{-1/2}σε^{-1/2}`），场按 `e^{−λ_i t}` 衰减，本征向量混合分量。

**配置**：

```
均匀 ε、σ（无界面）；σ = diag(σx,σy,σz) 或 R diag Rᵀ（绕 y 30°）
初始 E 沿主轴（对角情形）或沿本征向量（旋转情形）
探针：固定点，记录 Ex/Ey/Ez 时间序列
```

**分析**：

```
M1/M2：节点打印 vs Python expm（≤1e-10）
场：ln|E_comp| 斜率 vs −λ_i（≤2%）
旋转 σ：Ex/Ez 的耦合与本征向量方向
```

**已测基线**：M1 `4.7e-11`；M1_xz=0.04726（耦合生效）。

---

### L5 — 界面 + 张量 ε + σ

**物理**：s-pol 用 `ε_yy` 的 lossy Fresnel；板内 `δ_yy`；界面面电荷 `q_s=(ε2/σ2−ε1/σ1)J_n`（标量近似）或张量版的法向分量。

**配置**：

```
ε2 张量（对角/旋转）；sigma_function 只在材料 2；s-pol；16/32 ppw
```

**判据**：`R/T` ≤5%；`δ_yy` ≤2%；界面 Gauss 残差（见 L6）。

---

### L6 — 面电荷与连续性（ρ_medium）

**物理**（两导体弛豫，无需源）：

```
初始 t=0：两区域均匀 E0 ẑ（E_n≠0）
E1(t)=E0 e^{−t/τ1},  E2(t)=E0 e^{−t/τ2},  τ=ε/σ
q_s(t) = D2(t)−D1(t) = ε2 E0 e^{−t/τ2} − ε1 E0 e^{−t/τ1}
```

**配置**：

```
ε1,σ1（z<0）与 ε2,σ2（z>0）；界面 z=0
初始 E_z=E0（法向）；无源；记录界面节点 ρ_medium 与 E1/E2
```

**分析**：

```
L1: ρ(0) = q_s(0)/dz
L2: ρ_interface(t) vs q_s(t)/dz（≤2%）
L3: 界面 max|∇·D−(ρ_p+ρ_medium)| ~1e-13
```

**注意**：**不要**建"两导体 + 稳态 J_n"的算例——标准求解器没有 DC 电流源，E 会弛豫到零；稳态 q_s 需要先加外部电流源（功能开发）。

---

### L7 — 斜入射（p-pol，E_n≠0）

**物理**：p-pol 有 `E_n`，界面 `D_n` 跳变 = 波驱动的面电荷；Fresnel 用 `k_z`、`η` 的 p-pol 公式。

**配置**：

```
入射角 30°/45°；s-pol 对照；p-pol 用 E_n≠0
探针：反射/透射方向；DFT 分窗
```

**判据**：`R_p/T_p` ≤5%；界面 `ρ_medium` 与 `(D_n2−D_n1)` 的相位/幅度一致。

---

### L8 — 稳定性与鲁棒性

```
1. σ→∞（1e3 倍）：一步 E→0，无 NaN
2. 极端对比 ε2/ε1=100、σ2/σ1=1e6：长跑 1e4 步无 NaN、场有界
3. CFL 0.4 与 0.3：结果差 ≤1%（在收敛区）
4. no-op：σ=0 且无张量 → P0 逐位不变
5. 张量 σ 有负本征值（增益）：应警告/断言，不静默
```

---

### L9 — 粒子 + 介质（限制项）

当前 `ρ_medium` 与粒子总电荷**未接线**：Gauss 检查只含 ρ_medium。若跑"粒子在介质中"：

```
总自由电荷 = ρ_particle + ρ_medium
Gauss: ∇·D = ρ_p + ρ_medium
```

在接线前，这一项标注"未验证"。

---

### L10 — 明确不测（未交付/不支持）

- PML + 宏观介质（P8 不兼容）
- >2 种材料 / 多界面 / 三相交界
- AMR（lev>0）
- μ 3×3 张量（只有标量 μ）
- σ 的复 ε τ（Deinega，T2b）
- SIBC（明确拒绝：δ 被分辨，会绕过内部解）

---

## 3. 操作流程（每个算例照此执行）

### 3.1 构建

```bash
conda activate warpx-cpu-mpich-dev
cd /home/douliangduo/WarpX_fluid
cmake --build build_2d -j 8
```

### 3.2 Pre-flight（必做）

```bash
python3 benchmarks/cases/c03_medium/preflight.py benchmarks/cases/c03_medium/<case>/inputs_2d
```

检查五项：材料（σ 全域+界面 → FAIL）、路径衰减（>e^{-5} → FAIL）、三段时间、期望量级、探针 valid box。
**FAIL 的算例不进入 R/T 解读。**

### 3.3 运行

```bash
cd benchmarks/cases/c03_medium/<case>
OMP_NUM_THREADS=8 /home/douliangduo/WarpX_fluid/build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB \
    inputs_2d max_step=<N> <overrides>
```

- `dt` 从 `run.log` 读（矩阵验证必须用实际值）；
- 网格档位用 `amr.n_cell` 覆盖（480²/864²、16/32 ppw）；
- A/B 开关：`WARPX_YEE_CELL_EPS=1` → cell Kottke；不设 → node τ（默认）。

### 3.4 分析

```bash
python3 analyze.py <run_dir>        # R/T（DFT 分窗 + 通量）
python3 verify_divergence.py plt00010 --eps <ε> --rho-field rho_medium
```

- R/T 必须给出**窗口稳健性**（至少两个窗）；
- 能量账 `R+T+A`（A 解析）；
- Gauss 残差（内部/界面分别报）。

### 3.5 记录

每个算例记录：`R`、`T`、`A`、`δ_eff`、`v_φ`、`max|∇·D−ρ|`、`dt`、网格、窗、误差 vs 解析、是否 PASS。

---

## 4. 回归基线（当前默认：node τ）

| 算例 | 数值 | 备注 |
|---|---|---|
| P0（σ=0，512²） | `R_flux 0.0927 / T_flux 0.8906` | 逐位基线 |
| E3 480² | `R 0.141817 (−2.80%) / T 0.869037 (+1.75%)` | node τ |
| E3 864² | `R 0.144650 (−0.86%) / T 0.858631 (+0.53%)` | node τ，32 ppw |
| S5（有损板） | `T −0.215%` | 相位 −0.529° |
| E6（导体界面） | 深区 640² −1.58% / 1280² −0.61% | 量级二阶 |
| 界面 R/T（σ 只在材料 2） | `R 0.1552 (+4.9%) / T 0.870 (+2.1%)` | 共同偏置 ~5% |
| 张量+σ 内衰减 | `δ_eff=1.176 vs 1.187 (−0.91%)` | `v=c/√5` |
| 张量 σ 的 M1 | `4.7e-11`（vs Python） | bulk 旋转 |
| J0/J1 | `1e-10`（vs Python） | 张量 ε+σ |

---

## 5. 常见失败与定位

| 现象 | 先查 |
|---|---|
| `T≈0, R≈1` | σ 是否全域（preflight 材料项）；入射是否到达界面（路径预算） |
| `R_flux<0` | 真空侧入射/反射时间重叠；换反射窗/参考减 |
| NaN | σ 负本征值、ε 非 SPD、dt 不匹配 |
| Gauss 残差大 | ρ/⟨J_c⟩ 的 σ、ε、时间平均是否与 E 更新同源；σ 梯度处 |
| 张量结果错 | 节点 σ 的 off-diagonal、frame(n) 的切向/法向块、ε(I−C)ε⁻¹ 系数 |
| 与解析差 5% 常数偏置 | 归一化/窗口/参考跑；先看 R+T+A 是否同向偏 |

---

## 6. 参考

1. A. Taflove, S. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain Method*.
2. Meep, *Subpixel Smoothing*: https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/
3. Tidy3D, *Lossy Metal Subpixel Averaging*: https://tidy3d-help-center.simulation.cloud/flex_rf_gui/TCM/configuration/advanced/subpixel_averaging/LossyMetal.html
4. A. Deinega, I. Valuev, *Opt. Lett.* **32**, 3429 (2007): https://pubmed.ncbi.nlm.nih.gov/18059956/
5. A. Oskooi, A. Kottke, S. Johnson, *Opt. Lett.* **34**, 2778 (2009): https://math.mit.edu/~stevenj/papers/OskooiKo09.pdf
6. *Surface charge and surface current densities at material boundaries*, Am. J. Phys. **92**, 158 (2024).
7. *A continuity-preserving and divergence-cleaning algorithm ... in inhomogeneous media*, J. Comput. Phys. (2017).
8. NASA TM-2001-210862（有损介质 FDTD 算例）。
