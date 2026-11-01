# P0 之后的定量仿真测试清单（可复现）

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 日期：2026-11-01
> 范围：P0 = 各向同性介质、`maxwell_solver = yee`、`em_solver_medium = macroscopic`、subpixel 界面。
> 目的：在 P0 已通过"法向 Fresnel"之后，列出**还能定量测什么**、参考是什么（解析 / Meep / 自洽）、
> 验收多少、怎么复现。
> 原则：每个测试都给出「物理量 → 参考 → 验收 → 命令 → 期望值」，不写"跑通即可"。

---

## 0. 前提与已有结果

P0 已完成（commit `8851146` 修复后、`b2ccba0` 真空参考、`8d10a92` 加长域、`4ef439d` 文档）：

| 量 | 16 ppw（512²） | 32 ppw（1024²） | 理论 | 收敛 |
|---|---|---|---|---|
| `T_power` | 0.8917（+0.31%） | 0.8851（−0.43%） | 8/9 = 0.8889 | 稳 |
| `R + T` | 1.0223（+2.2%） | 1.0007（+0.07%） | 1 | O(h²) |
| `R_flux` | 0.1335（+20.1%） | 0.1171（+5.40%） | 1/9 = 0.1111 | ≈O(h²) |
| 点值 `R` | 0.1307（+17.6%） | 0.1155（+3.99%） | 1/9 | ≈O(h²) |
| 相位（z=1.0） | — | −14.2° | 0° | 色散量级 |

**已确立的口径**：`T_power` ≥16 ppw ≤1%；`R+T` ≥32 ppw ≤1%；`R_flux` 按 O(h²) 收敛。
**已知限制**：PML 与 macroscopic 不兼容（P8）→ 用加长域；界面逐点场至多一阶（Meep limitation 4）。

---

## 1. 测试矩阵（按成本排序）

| # | 测试 | 验证的物理 | 参考 | 验收 | 成本 |
|---|---|---|---|---|---|
| T1 | 法向 Fresnel R/T | 界面阻抗、能量守恒 | 解析 | 已通过（见 §0） | 已完成 |
| T2 | 斜入射 s/p Fresnel + Brewster + 临界角 | 角度依赖、偏振 | 解析 | R/T ≤1%（32 ppw）；Brewster 角 ±0.5° | 中 |
| T3 | Snell 定律（折射角） | k 矢量匹配 | 解析 | 折射角 ≤0.5° | 低 |
| T4 | 薄板 Fabry–Perot | 多次反射、相位 | Airy 公式 | 峰位 ≤1%、峰高 ≤2% | 中 |
| T5 | 均匀介质数值色散 ω(k) | Yee 色散 | 解析（Yee） | 相速度 ≤0.5% | 低 |
| T6 | 能量守恒 / 损耗衰减 | 无源/有损 | 解析 | 无损 R+T=1（32 ppw ≤1%） | 低 |
| T7 | Meep ring resonator 频率 | subpixel 的二阶性 | Meep | 频率误差 ≤0.5%（同 ppw） | 中 |
| T8 | Meep waveguide bend 通量 | 曲面界面、散射 | Meep | 通量 ≤2% | 中 |
| T9 | 2D 介质圆柱 Mie / RCS | 曲面、远场 | 解析 Mie | RCS ≤2%（同 ppw） | 中 |
| T10 | 退化/等价性回归 | 代码路径 | 自洽 | 逐位（ε1=ε2、tensor=scalar） | 低 |

---

## 2. 各测试的详细方案

### T1 法向 Fresnel（已完成，作为基线）

- 算例：`ε1=1`、`ε2=4`、界面 z=0、`maxwell_solver=yee`、`macroscopic_subpixel=1`；
- 方法：**真空参考跑**（同源、同网格）+ 双探针分频域相减；坡印廷通量 `S_z=(Ex·By−Ey·Bx)/μ0`（各向同性时退化为 `−Ey·Bx/μ0`）；
- 理论：`R=1/9`、`T_power=8/9`；
- 命令：见 §3；期望：§0 表。
- **注意**：PML 不可用（P8），必须用加长域 `[-8,8]`，窗口避开 PEC 反射。

### T2 斜入射 s/p Fresnel + Brewster + 临界角

- **物理**：s/p 反射系数随角度的变化；Brewster 角（`tanθ_B=n2/n1`，n2/n1=2 → 63.4°）；从密到疏（ε1=4→ε2=1）的临界角（`sinθ_c=1/2` → 30°）与全内反射。
- **方法**（二选一）：
  1. **周期边界 + 平面波源**：`boundary.field_lo/hi = periodic`（y 方向周期），入射波用 `k=(kx,0,kz)` 的解析式初始化；DFT 取 R/T；
  2. **真空参考 + 双探针**（P0 方法的斜入射版）：入射侧两个探针做行波分解，`k` 用**真空跑测得的 k_num**（不用 `ω/c`）。
- **验收**：32 ppw 下 `|R_s|,|R_p| ≤1%`；Brewster 角位置 ±0.5°；临界角 ±0.5°；全内反射时 `|R|≈1`（≤0.5%）。
- **成本**：中；需要斜入射初始化与斜探针。

### T3 Snell 定律

- **物理**：`n1 sinθ1 = n2 sinθ2`；用空间 FFT（plotfile）测两侧的 `k` 方向。
- **方法**：两个 plotfile（入射侧、透射侧），对 `Ey(x,z)` 做 2D FFT，取峰值 `k` 矢量，算角度。
- **验收**：折射角 ≤0.5°（32 ppw）。
- **成本**：低；需要 plotfile 输出（`diag1.fields_to_plot = Ey`）。

### T4 薄板 Fabry–Perot

- **物理**：厚度 d 的板（ε2）两侧为 ε1，透射谱的 Airy 振荡：`T = 1/(1+F sin²(δ/2))`，`δ=2k2d`，`F=4R/(1−R)²`。
- **方法**：宽带脉冲 + 透射侧 DFT，扫 d（或扫频）；比较峰位/峰高。
- **验收**：峰位 ≤1%、峰高 ≤2%（32 ppw）。
- **成本**：中；与 T2 共享斜入射/DFT 工具。

### T5 均匀介质数值色散

- **物理**：Yee 的 `ω(k)`；相速度/群速度。
- **方法**：均匀 ε2 介质（无界面）+ 宽脉冲；探针测群速度；或空间 FFT 测 `k(ω)`；与 Yee 解析色散比。
- **验收**：相速度 ≤0.5%（~30 ppw）。
- **成本**：低；已有 `verify_medium_speed.py` 可复用。

### T6 能量守恒 / 损耗

- **无损**：`R+T=1`（32 ppw ≤1%）；
- **有损**（P2 之后）：均匀导体 `ρ ∝ e^{−t/(ε/σ)}`，时间常数 ≤2%；能量单调衰减；
- **成本**：低。

### T7 Meep ring resonator 频率收敛

- **物理**：subpixel 的二阶性（Meep 的经典测试：ring resonator 的共振频率随分辨率二阶收敛）。
- **方法**：
  - Meep：2D ring（`ε=12` 硅环 + 空气），偶极源，`harminv` 提取频率，分辨率 10/20/40（ppw）；
  - WarpX：同一几何（level set/`interface_function` 圆环）、同源、同 ppw，DFT/峰值提取频率；
  - 比较两者对高分辨率参考的误差。
- **验收**：同 ppw 下频率相对误差 ≤0.5%；两者都呈二阶收敛。
- **成本**：中；需要 Meep 脚本 + WarpX 圆环算例。

### T8 Meep waveguide bend 通量

- **物理**：曲面界面 + 散射；90° 弯波导的透射率。
- **方法**：Meep 教程的 bend 算例；WarpX 同几何；比较输出端通量（`S_z` 积分）。
- **验收**：通量 ≤2%（同 ppw）。
- **成本**：中。

### T9 2D 介质圆柱 Mie / RCS

- **物理**：曲面界面、远场散射。
- **方法**：复用 `mie_cylinder.py`（把 PEC 换成 ε2=4 的介质，解析 Mie 用贝塞尔函数）；`analyze_mie_convergence.py` 做收敛；与 Meep 同 ppw 对照。
- **验收**：RCS ≤2%（同 ppw）；近场 L2 二阶收敛。
- **成本**：中；已有 PEC 版工具，介质版需改解析解。

### T10 退化/等价性回归

- `ε1=ε2` 倾斜界面：tensor 与 scalar 逐位（`max|diff|=0`）；
- `epsilon = ε` 与 `epsilon_tensor_xx=yy=zz=ε`：逐位；
- 各向同性 Fresnel：`R_flux=0.1335`、`T_flux=0.8513` 不变；
- **成本**：低；应进 CI。

---

## 3. 复现协议

### 3.1 环境

- WarpX：本仓库 `WarpX_fluid`，commit 记录在报告头部（P0 基线 `8851146` + 后续）；
- 构建：`cmake -S . -B build -DWarpX_DIMS=2 -DWarpX_PYTHON=ON`；`cmake --build build -j`；
- Python：`numpy`、`matplotlib`（分析脚本）；
- Meep：版本以 `meep --version` 为准（本机文档为 1.17.1）；Meep 脚本放 `benchmarks/cases/c02_0_pec/meep/`。

### 3.2 命令模板

```bash
# T1 法向 Fresnel（真空参考法）
cd benchmarks/cases/c02_0_pec
python3 verify_fresnel_rt.py <run_root> --eps2 4.0 --ncell 512 --nsteps 1400 \
    --tmax-inc 2.0e-8 --tmax-trans 2.6e-8 --make --execute
# 期望：R_flux≈0.1335, T_flux≈0.8513（16 ppw）；R+T≈0.9847
#       32 ppw：R_flux≈0.1171, T_flux≈0.8778, R+T≈0.9949
```

```bash
# T2/T3 斜入射（T2 用周期边界 + DFT；T3 用 plotfile + 2D FFT）
# T4 薄板：同 T2 的 DFT 工具，扫厚度
# T5 色散：verify_medium_speed.py（宽脉冲 + 包络峰）
```

### 3.3 期望值与容差

- 每个测试脚本输出：`R/T`、`R+T`、相位、收敛率；容差见 §1 表；
- 回归：T1 的数字必须与 §0 表一致（逐位/≤1e-14）。

### 3.4 数据与脚本位置

```
benchmarks/cases/c02_0_pec/
  verify_fresnel_rt.py      # T1（真空参考 + 坡印廷通量）
  analyze_oblique.py        # T2 参考
  verify_medium_speed.py    # T5
  mie_cylinder.py           # T9 解析
  analyze_mie_convergence.py# T9 收敛
  tau_transform_test.py     # τ 单元测试（P1）
  generalized_fresnel.py    # 张量参考（P1-d）
```

---

## 4. 优先级与建议顺序

```
T10（回归，分钟级） → T5（色散，低） → T1（基线，已有） →
T2（斜入射 Brewster/全反射，中） → T3（Snell） → T4（Fabry–Perot） →
T7（Meep ring） → T8（Meep bend） → T9（介质 Mie/RCS）
```

理由：先自洽/解析（便宜、能立刻排除代码问题），再上 Meep（跨代码对照成本高）；T2 是斜入射的入口，T3/T4 复用它的工具。

---

## 5. 限制与注意事项

1. **PML 不可用**（P8）：`em_solver_medium=macroscopic` + PML 会 abort；用加长域 + 窗口；
2. **界面逐点一阶**（Meep limitation 4）：验收用**界面外/积分量**（R_flux、通量、RCS、共振频率），不用界面点值；
3. **R 是 O(h²)**：16 ppw 的 +20% 是正常离散误差；验收按收敛率，不按 16 ppw 的 1%；
4. **subpixel 仅适用无损、无频散介质**（Meep limitation 1）：有损/色散（P2）另算；
5. **Meep 对照要同 ppw、同 subpixel 开关、同源**；Meep 的 Yee 时间步/CFL 与 WarpX 不同，比较的是**收敛值/频率/通量**，不是逐点场；
6. 所有期望值以"高分辨率参考"或解析值为准，脚本里硬编码期望值并做断言。

---

## 6. 参考文献

1. A. Farjadpour et al., *Opt. Lett.* 31, 2972 (2006)（subpixel）。
2. C. Kottke, A. Farjadpour, S. G. Johnson, *Phys. Rev. E* 77, 036611 (2008)，DOI 10.1103/PhysRevE.77.036611。
3. Meep 文档：Subpixel_Smoothing（limitation 1/4）、ring resonator / waveguide bend 教程。
4. Meep 源码：`anisotropic_averaging.cpp`（`eff_chi1inv_row`）。
5. 本项目：`技术路径.md` §7（验证与基准）、`verify_fresnel_rt.py`。

*本清单为 P0 之后的定量测试方案；T1 已有数据，T2–T9 为待做项，按 §4 顺序推进。*
