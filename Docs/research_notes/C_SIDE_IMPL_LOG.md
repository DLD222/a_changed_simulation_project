# C 侧实施记录（S0'–S3' 完成）

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 日期：2026-10-31 · 分支 `feat/ect-conformal-E`
> 依据：`C_SIDE_DESIGN2.md`（实现）、`C_SIDE_RESEARCH2.md`（开工条件）、`Docs/C_SIDE_REVIEW.md`（审查）

---

## 1. 已落地的代码

| 文件 | 改动 |
|---|---|
| `Source/Fields.H` | 新增 `FieldType::dual_A_y` / `dual_L_x` / `dual_L_z` |
| `Source/WarpX.H` | 新增 `static inline int ect_conformal_E = 0` |
| `Source/WarpX.cpp` | 解析 `algo.ect_conformal_E`（0/1/2，含 2D-only 与 ECT-only 断言）；ECT 分配段按 `Ey/Bx/Bz` 交错分配三套对偶量 |
| `Source/EmbeddedBoundary/EmbeddedBoundaryInit.{H,cpp}` | 新增 `ComputeDualGeometry`（`constant_geometry` 与真实几何两条分支，方案 B） |
| `Source/Initialization/WarpXInitData.cpp` | `InitializeEBGridData` 的 ECT 分支调用，并在 `verbose` 时打印 min/max |
| `Source/FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H` | 声明 `EvolveECartesianECT` |
| 同上 `EvolveE.cpp` | 分派（`ECT && ect_conformal_E>0`）+ 内核实现（三分量，2D XZ） |
| `benchmarks/cases/c02_0_pec/verify_l2a.py` | 判据加入"噪声地板"（零分量用其它分量的 1e-13 作绝对容差） |

开关语义：`algo.ect_conformal_E = 0`（关，旧路径）/ `1`（常数几何，调试）/ `2`（真实几何）。

内核（2D XZ）：

```
Ex(i,j) += −c²dt [By(i,j) − By(i,j-1)] / Lz(i,j)
Ez(i,j) += +c²dt [By(i,j) − By(i-1,j)] / Lx(i,j)
Ey(i,j) += +c²dt [Bx(i,j)Lx(i,j) − Bx(i,j-1)Lx(i,j-1)
                  −Bz(i,j)Lz(i,j) + Bz(i-1,j)Lz(i-1,j)] / Ay(i,j)
```

几何量（方案 B，`v = volfrac`）：

```
Ay(i,j) = ¼ dx dz [v(i-1,j-1)+v(i,j-1)+v(i-1,j)+v(i,j)]
Lx(i,j) = ½ dx [v(i-1,j) + v(i,j)]        # 与 Bx 同位，IntVect (1,0)
Lz(i,j) = ½ dz [v(i,j-1) + v(i,j)]        # 与 Bz 同位，IntVect (0,1)
```

---

## 2. 验证结果

### 2.1 L2a（规则域，EB 在域外）

| 配置 | Ey rel | Bx rel | 能量 rel | 判定 |
|---|---|---|---|---|
| 常数几何（`=1`） | 1.6e-15 | 5.2e-16 | 2.6e-15 | **PASS** |
| 真实几何（`=2`） | 1.9e-15 | 6.2e-16 | 2.6e-15 | **PASS** |

几何量打印（64², dx=dz=0.25）：`A_y = 0.0625..0.0625`（= dx·dz）、`L_x = L_z = 0.25`（= dx = dz），
真实几何分支在规则域给出**精确常数**。

### 2.2 重要修正：不可能"逐位"，只能到机器精度

设计文档里的 C-1 判据写的是 `max|d| = 0`。实际实现后 `max|d| ≈ 1e-9`（相对 1.6e-15），
原因不是索引错，而是**浮点运算顺序不同**：

```
Yee    : (Bx(i,j) − Bx(i,j-1)) / dz
共形   : (Bx(i,j)·dx − Bx(i,j-1)·dx) / (dx·dz)
```

两者数学等价，但先乘后减再除 ≠ 先减后除。**任何"乘开放长度再除对偶面积"的共形形式都不可能
与差商形式逐位相同**。因此 C-1 的正确表述是"到机器精度（`rel ≤ 1e-14`）"，本文与后续文档统一采用。

副产物：`Bz` 在 TE 算例中峰值只有 `1e-19`（纯舍入噪声），旧的相对判据会给出 `rel = 1.0` 的假 FAIL；
`verify_l2a.py` 已加噪声地板修正。

### 2.3 L3 45° TE（真实几何）

| 判据 | 旧（掩膜 Yee） | 新（共形 E） | 变化 |
|---|---|---|---|
| 平均 `\|r+1\|` | 0.0427 | 0.0459 | 略差 |
| 斜率误差（角度） | −2.39% | **−1.51%** | **改善** |
| 单点起伏 | −0.991 .. −1.169 | −0.996 .. −1.174 | 持平 |
| 交叉偏振 | 0 | 0 | — |

结论：在 45° 斜面上，共形 E 的"去掉阶梯误差"与"一阶几何引入的新误差"量级相同（RESEARCH2 Q12 的预测），
**净效果不显著**。要看到实质提升需要方案 C（level-set 二阶几何）。

---

## 3. 关键发现：米氏"爆炸"与 C 侧无关

原本的判断（`HANDOFF.md` §3）是"曲边界极小胞元上 B 共形/E 阶梯不一致 → 必须做 C 侧"。
实测推翻了它。全部在 80²、域 `[-4,4]²`、400 步下：

| # | 配置 | step 200 的 `max\|E\|` | 结论 |
|---|---|---|---|
| 1 | ECT + EB 圆柱 + 掩膜 Yee（原始） | ~1e306（数十步即爆） | 爆 |
| 2 | ECT + EB 圆柱 + 共形 E（常数几何） | 1.03e308 | 爆 |
| 3 | ECT + EB 圆柱 + 共形 E（真实几何） | step 64 已 5.8e307 | 爆（更早） |
| 4 | ECT + **EB 移出域外** + 共形 E | 1.36e308 | **爆** |
| 5 | **Yee + EB 移出域外 + PML** | 1.67e308 | **爆** |
| 6 | **Yee + EB 移出域外 + PEC** | 1.46e308 | **爆** |
| 7 | 同 4，但**源幅值置零** | 2.88e6（400 步不变） | **稳定** |

第 5/6 组证明：**不用 EB、不用 ECT、只用 Yee，同样爆炸**。第 7 组证明：源一旦置零就完全稳定。

**真正根因（同一轮内已定位）**：`amr.max_grid_size = 64` 而 `n_cell = 80` ⇒ 域被切成 4 个 grid。
CW 激光源**跨 grid 边界时不稳**。判据是单/多 grid 的分界，与 EB、ECT、C 侧、PML 全部无关：

| 配置（80²，无 EB，400 步） | `max_grid_size` | 末帧 `max\|E\|` |
|---|---|---|
| ECT，conformal 关 | 64（4 grid） | **1.76e308** |
| ECT，conformal 开 | 64（4 grid） | 1.36e308 |
| ECT，conformal 关 | **80（1 grid）** | **1.25e6** |
| ECT，conformal 开 | **80（1 grid）** | **1.25e6** |
| Yee | **80（1 grid）** | **1.25e6** |
| ECT + **EB 圆柱**，conformal 关 | **80（1 grid）** | **1.92e6** |
| ECT + **EB 圆柱**，conformal 开 | **80（1 grid）** | **1.86e6** |

⇒ 这同时解释了历史现象：40²/64² 时 `max_grid_size=64 ≥ n_cell`（单 grid）→ 稳定；
80²/160² 时被切块 → 爆炸。**与"极小 volfrac 胞元"无关**。

修复：`make_mie_cases.py` 的模板改为 `amr.max_grid_size = {n}`（单 grid）。

> 遗留问题（独立 bug，不在 C 侧范围）：CW 激光源 + 多 grid 的不稳定机理未查，
> 生产算例若必须分块需先解决（或用脉冲/初值替代 CW 源）。

### 3.1 米氏收敛阶基线（单 grid + 后处理修正后）

后处理另有两个 bug，本轮一并修掉（`analyze_mie_convergence.py`）：

1. **相位参考错**：DFT 用 `exp(i·ω·j·T/8)`，但帧间隔是整数个时间步（40² 时 `interval = 2`，
   与 T/8 差 12%），累积相位误差 ~34°。改为读 plotfile Header 第 5 行的物理时间。
2. **未归一化**：数值场是 V/m（`E0 = 1e6`），解析场幅值为 1，L2 误差被 `E0` 支配（恒为 1e6）。
   改为除以 `my_constants.E0`。

修正后（TE，近场环 `1.5R–2.5R`，8 帧 DFT）：

| 配置 | 近场 L2 @80² | @160² | 近场阶数 | RCS 阶数 |
|---|---|---|---|---|
| 掩膜 Yee（`ect_conformal_E = 0`） | 0.794 | 0.327 | **1.28** | **1.76** |
| C 侧（`= 2`，方案 B 一阶几何） | 0.932 | 0.440 | 1.06 | 1.77 |

结论：**一阶对偶几何与阶梯同阶**（近场 ~1 阶），C 侧要在收敛阶上兑现必须换方案 C
（level-set 精确几何，DESIGN2 §1.5）。RCS 两者都接近 2 阶，因为远场外推平滑了边界误差。
这两个数字就是 S10' 的对照基线。

**对项目的影响**：
1. `HANDOFF.md` §3 的"正解 = C 侧"结论作废，需改写；
2. 米氏收敛阶测试**当前不可用**，但不是被 C 侧阻塞；
3. C 侧的价值回到"精度与守恒"（TM 能量增长、界面附近误差），而不是"修稳定性"。

---

## 4. 下一步（按优先级）

1. **定位激光源的不稳定**（阻塞米氏）：
   - 用 `sin(omega*t)` 替代 `cos`（去 DC），加平滑开启 `*(1-exp(-t/tau))`；
   - 对比 40²/64²/80²/128² 找触发阈值（同域、只改 `n_cell`）；
   - 检查 `FieldFunctionLaserProfile` 是否每步 `+=` 注入（若是，源平面内部会出现反馈式叠加）。
2. **方案 C（level-set 二阶对偶几何）**：45° 上把 `0.0459` 压下去，这是 C 侧真正兑现精度的地方。
3. **对偶稳定化**（`A_stab`）：本轮未做；在真实几何上实测 `A_y` 最小到 0、`L` 最小到 0，
   covered 点靠内核的 `<= 0` 分支置零，cut 点的小值目前**无保护**。
4. **TM 长跑 1e4 步**：验证能量漂移是否从 9.5e-3 降到 ≤1e-4（C 侧的直接目标之一）。
5. **3D**：内核目前 2D-only（`#if defined(WARPX_DIM_XZ)`），3D 构建会在 `WarpX.cpp` 断言处 abort。

---

## 5. 复现命令

```bash
cd <repo>
cmake --build build_2d -j 8

SRC=benchmarks/cases/c02_0_pec
BIN=build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB

# L2a：常数几何 / 真实几何，均应与 Yee 到机器精度
mkdir -p /tmp/c1/ect_real && cd /tmp/c1/ect_real
sed 's/^\(algo\.maxwell_solver *= *ect\)$/\1\nalgo.ect_conformal_E = 2/' \
    <repo>/$SRC/inputs_ect_2d_noeb > inputs
OMP_NUM_THREADS=8 <repo>/$BIN inputs > run.log 2>&1
python3 <repo>/$SRC/verify_l2a.py . /tmp/l2/yee_noeb

# L3 45° TE
mkdir -p /tmp/c1/l3_45te && cd /tmp/c1/l3_45te
sed 's/^\(algo\.maxwell_solver *= *ect\)$/\1\nalgo.ect_conformal_E = 2/' \
    <repo>/$SRC/inputs_2d_oblique_te_45deg > inputs
OMP_NUM_THREADS=8 <repo>/$BIN inputs > run.log 2>&1
python3 <repo>/$SRC/analyze_oblique.py . 45 te
```
