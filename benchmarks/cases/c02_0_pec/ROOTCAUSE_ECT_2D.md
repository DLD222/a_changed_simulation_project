# ECT 场传播根因报告（2D XZ）

> 算例：`benchmarks/cases/c02_0_pec/`
> 代码：`<repo>`（WarpX 26.10 分支 + `Source/MomentsSolver`）
> 对照基线：`~/WarpX_changed`（2026-10-31 全新 clone）
> 日期：2026-10-31

---

## 0. 结论摘要（先看这四条）

1. **fork 没有回归。** `Source/` 与上游 clone 逐字节一致（唯一差异是新增的 `Source/MomentsSolver/`）。
   `diff -rq` 全仓库仅 1 处 Source 差异（新增目录），无任何 WarpX 原有文件被修改。
2. **fork 的 2D ECT 数值正确。** 上游回归算例 `embedded_boundary_rotated_cube`（2D XZ + ECT）
   在本 fork 上复算，场 checksum 与上游基准 JSON 匹配到 ~1e-9 相对误差（见 §3.2）。
3. **根因是上游 2D ECT 的实现边界**：`EvolveRhoCartesianECT` 在 `WARPX_DIM_XZ` 下只计算 `Rhoy`
   （TM 模：Ex/Ez → By），**不计算 `Rhox` 与 `Rhoz`**，因此 **TE 模（Ey → Bx/Bz）在 2D XZ 下不被支持**。
   这不是 bug，而是当前上游 2D 实现的适用范围（见 §4）。
4. **换到 TM 模后一切正常**：同一斜界面/正入射算例，TM 偏振下 PEC 反射系数
   **r = −0.9993**（期望 −1，误差 0.07%），入射/反射时刻与解析值吻合（见 §3.4）。

**一句话**：问题不在 fork，而在于最初的 C2-0 用了 **TE 偏振（Ey 出平面）**，落在上游 2D ECT 的
未实现分支上。2D XZ ECT 只做 TM。要测 TE，必须用 3D，或扩展 `Rhox`/`Rhoz` 的 XZ 实现。

---

## 1. 问题陈述（复现到的两个异常）

在 `c02_0_pec` 的 C2-0 基线（倾斜 PEC 平面 + 平面波）上，最初观测到：

| 现象 | 配置 | 观测 |
|---|---|---|
| A | 2D XZ + ECT + **激光**（polarization `0 1 0`，即 Ey） | 探针处 Ey 始终为 0 |
| B | 2D XZ + ECT + **场初始化**（Ey + Bx） | Ey 单调上升、59 ns 仍未过峰（应为 10 ns 过峰） |
| C | 2D XZ + **yee** + 同源 | Ey ≈ 9.1e5 V/m，正常传播 |

A/B 只在 ECT 下出现；yee 正常。由此把问题锁定在 ECT 路径。

---

## 2. 证据链总览

```
现象 A/B（ECT 下零场/异常传播）
        │
        ├── 排除 fork 回归 ──► §3.1 源码逐字节一致
        │                    └► §3.2 rotated_cube checksum 对齐上游
        │
        └── 定位到 2D XZ ECT 的 Rho 分支 ──► §4.1 只算 Rhoy
                                    ├── §4.2 TE 需要 Rhox/Rhoz（未算）
                                    └── §3.3/§3.4 TE 失败 vs TM 成功
```

---

## 3. 证据

### 3.1 fork 未修改任何 WarpX 原源码

命令与结果：

```text
$ diff -rq ~/WarpX_changed/Source ~/WarpX_fluid/Source
Only in <repo>/Source: MomentsSolver
（共 1 行输出）

$ diff -rq ~/WarpX_changed ~/WarpX_fluid \
      -x build -x build_2d -x .git -x MomentsSolver -x benchmarks -x diags -x Backtrace
（仅剩运行时产物：warpx_used_inputs / Backtrace.0.0 / 既有软链接）
```

结论：`Source/FieldSolver/**`（含 `EvolveB.cpp` 的 ECT、`EvolveECTRho.cpp`、`EvolveE.cpp`）
与上游 clone **完全一致**。因此 fork 不可能在 ECT 代码里引入回归——行为差异只能来自
“构建/输入/上游实现边界”。

### 3.2 rotated_cube checksum 与上游基准一致

上游回归算例：`Examples/Tests/embedded_boundary_rotated_cube/inputs_test_2d_embedded_boundary_rotated_cube`
（2D XZ，`algo.maxwell_solver = ect`，PEC 旋转立方体，By 余弦模初始化）。

运行后 plotfile：`diags/diag1000068`。用自写的最小 plotfile 读取器
（`benchmarks/plotfile_checksum.py`，解析 `Header` + `Level_0/Cell_D_00000`，
按分量求 `Σ|Q|`）得到：

| 分量 | 本 fork 复算 | 上游基准 `test_2d_embedded_boundary_rotated_cube.json` | 相对差 |
|---|---|---|---|
| Ex | 7.8275781e+04 | 7.8275781e+04 | ~1e-9 |
| Ey | 0 | 0 | — |
| Ez | 3.1681636e+04 | 3.1681636e+04 | ~1e-9 |
| Bx | 0 | 0 | — |
| By | 2.7221007e-04 | 2.7221007e-04 | ~1e-9 |
| Bz | 0 | 0 | — |

结论：**fork 的 2D ECT 与上游基准数值一致**（差异仅在浮点求和顺序/AMReX 版本量级）。
特别地，E 场分量只有 **Ex、Ez** 非零、B 场只有 **By** 非零——这正是 **TM 模**。

### 3.3 关键代码事实：`EvolveRhoCartesianECT` 在 XZ 下只算 `Rhoy`

文件：`Source/FieldSolver/FiniteDifferenceSolver/EvolveECTRho.cpp`（`EvolveRhoCartesianECT`，约 111–146 行）：

```cpp
amrex::ParallelFor(trhox, trhoy, trhoz,
  // Rhox (Bx 面):
  [=]{ if (Sx<=0) return;
#ifndef WARPX_DIM_XZ
       Rhox = (Ey*ly - Ey*ly + Ez*lz - Ez*lz) / Sx;     // XZ 下整段不执行
#endif
  },
  // Rhoy (By 面):
  [=]{ if (Sy<=0) return;
#ifdef WARPX_DIM_XZ
       Rhoy = (Ez*lz - Ez*lz + Ex*lx - Ex*lx) / Sy;      // XZ 分支
#elif defined(WARPX_DIM_3D)
       Rhoy = (...);
#endif
  },
  // Rhoz (Bz 面):
  [=]{ if (Sz<=0) return;
#ifndef WARPX_DIM_XZ
       Rhoz = (Ex*lx - Ex*lx + Ey*ly - Ey*ly) / Sz;     // XZ 下整段不执行
#endif
  });
```

（为便于阅读省略了 stencil 索引，语义与原码一致。）

要点：

- `Rhox` 与 `Rhoz` 的赋值被 `#ifndef WARPX_DIM_XZ` 包住 ⇒ **2D XZ 下不写这两个数组**。
- 只有 `Rhoy` 在 `WARPX_DIM_XZ` 分支里被计算。
- 而 `EvolveBCartesianECT`（`EvolveB.cpp`）对 `idim = 0,1,2`（Bx/By/Bz）统一使用
  `Rho[idim]` 做 `B -= dt * rho_enl`。于是 **Bx、Bz 的 Rho 是未更新的旧值/初值**。

### 3.4 TE 失败 vs TM 成功的对照实验

全部在 `benchmarks/cases/c02_0_pec/` 下，2D XZ、`algo.maxwell_solver = ect`、
PEC 平面（`warpx.eb_implicit_function`）、GPU-free（CPU/OMP）、λ = 1 m、域 30 m、
探针 `(x,z) = (0,-5)`（距初始脉冲中心 3 m，距平面约 5 m）。

**(a) TE 偏振（Ey 出平面，Bx 入平面）——失败**

```text
inputs_2d      : Ey = E0·exp(-((z-zc)/w)^2),  Bx = -Ey/c
观测            : Ey(t) 从 1.8e4 单调升到 8.8e5（59 ns），不过峰、无反射
预期            : 10 ns 过峰（3 m / c），43 ns 反射峰
```

解释：TE 需要 `∂Ey/∂t = c²(∂Bx/∂z − ∂Bz/∂x)` 与 `∂Bx/∂t = ∂Ey/∂z` 闭合。
但 `Bx` 的 Rho（`Rhox`）在 XZ 下不更新 ⇒ Bx 冻结在初值。
于是 `∂Ey/∂t = c²∂Bx/∂z`（Bx 空间梯度固定不变）变成**恒定源**，Ey 只增不减——正是观测到的单调上升。

**(b) TM 偏振（Ex/Ez 入平面，By 出平面）——成功**

```text
inputs_2d_tm   : Ex = E0·exp(-((z-zc)/w)^2),  By = Ex/c
斜界面 α=30°   : max|Ex| = 9.995e5；入射 10.02 ns 到探针；产生 Ez（斜反射）
```

**(c) TM + 正入射（α=π/2，平面 z=0）——定量 r = −1**

```text
inputs_2d_tm_flat : alpha = pi/2 → eb = z，平面 z=0，法向 z
E0 = 1.0e6 V/m

入射峰 : t = 10.02 ns ,  Ex = +1.0000e6   （3 m / c = 10.0 ns ✓）
反射峰 : t = 43.49 ns ,  Ex = -9.9933e5   （10 + 2×5 m/c = 43.3 ns ✓）
r = 反射/入射 = -0.9993                     （PEC 期望 -1，误差 0.07%）
```

结论：TM 偏振下 ECT 的场传播与 PEC 反射**完全正确**（r 误差 0.07%，时刻误差 < 0.2 ns）。
这同时反证：ECT 机器本身没问题，问题只出在 TE 分支。

**(d) ECT + 激光（现象 A）**

激光 `polarization = 0 1 0`（即 Ey，TE）注入。`LaserParticleContainer::Evolve` 正常运行
（profiler 可见），但由于 TE 的 Bx/Bz 不更新，注入的场无法形成自洽的电磁波，
叠加激光注入机制的边界处理，探针处表现为 0。改用 yee（`inputs_yee`）则 Ey ≈ 9.1e5 正常。
这不是激光本身的 bug，而是同一 TE 限制的另一种表现。

---

## 4. 根因分析

### 4.1 为什么 2D XZ 只做 TM

2D XZ（`WARPX_DIM_XZ`）中，坐标是 x、z，y 为出平面方向。场分量 Ex、Ez 在平面内，
Ey 出平面；Bx、Bz 在平面内，By 出平面。一个沿 z 传播的平面波可分解为：

- **TM 模**：E 在平面内（Ex），B 出平面（By）。Faraday 只需 `Rhoy = ∮E·dl / Sy`，
  而 `Sy`（y 面 = 元胞本身）在 2D 有明确的面积分数。
- **TE 模**：E 出平面（Ey），B 在平面内（Bx、Bz）。需要 `Rhox`、`Rhoz`——
  它们对应 x 面 / z 面，而在 2D 中这些“面”退化为线段，其 cut-cell 面积/环量语义
  与 3D 不同，上游尚未实现（`#ifndef WARPX_DIM_XZ`）。

上游的 2D ECT 回归算例（rotated_cube、以及 `embedded_boundary_diffraction`）**全部是 TM**，
与“只做 TM”一致。这不是疏漏被发现，而是当前实现的既定适用范围。

### 4.2 相关实现约束（同源）与评审 D2 的呼应

`EvolveBCartesianECT` 仅支持 `WARPX_DIM_3D`/`WARPX_DIM_XZ`；ECT 几何量只在 `maxLevel()` 分配
（`WarpXInitData.cpp`）；E 更新（`EvolveE.cpp`）在 ECT 下与 Yee **同码**、仅掩膜跳点，无共形 E 更新。
评审 D2 指出 E 场共形 Ampere 需要**对偶面/对偶边**（现成 `face_areas`/`edge_lengths` 是 B 面/E 边，
不是 E 的对偶量）。本报告从另一侧印证：即便只为让 **Bx/Bz** 在 2D 正确演化，也要为 x/z 面补出
2D 语义的 `Rhox`/`Rhoz`（= 对偶几何 + 环量），工作量不小。

---

## 5. 影响与对策

### 5.1 对当前验证计划的影响

- **C2-0 基线**：改用 TM（Ex/By）即可跑通，已得 r = −0.9993，可作为 PEC cut-cell 几何自检基线。
- **TE（Ey 出平面）算例**：2D XZ ECT 下无效，必须用 3D（`WARPX_DIM_3D`，3D3V），或等 2D Rhox/Rhoz 实现。
- **激光 + ECT**：需用 TM 偏振（Ex/Ez），或转 3D。

### 5.2 三条可选路线

| 路线 | 内容 | 成本 | 风险 |
|---|---|---|---|
| R-A | 验证全部改用 **TM**（2D）或 **3D** | 低 | 3D 计算量大；TM 覆盖不了 TE 类问题 |
| R-B | 在 `EvolveRhoCartesianECT` 补 **2D 的 Rhox/Rhoz** | 中 | 需要正确 2D 对偶环量语义；必须与 3D 分支一致性测试 |
| R-C | 直接上 **E 场共形 Ampere + 对偶几何**（评审 D2） | 高 | 三件新工作；但这是完整 cut-cell Maxwell 的必经之路 |

**建议**：短期用 R-A 建立 TM 基线并推进介质/流体验证；中期做 R-B（补 2D TE）以支持
任意偏振的 2D 算例；R-C 作为 PEC→介质完整共形的核心实现，单独排期。

### 5.3 对 10 算例计划的具体修订

- C2-0：偏振改 TM（或 3D 跑 TE）。C2（介质 cut-cell，M2）：2D 限 TM，否则 3D。
- C4（Mie）：2D 圆柱建议 TM（按 Bohren–Huffman 约定核对），或直接用 3D 球（ECT 3D 无此限制）。

---

## 6. 复现步骤

```bash
# 0) 构建（2D, EB, ECT）
cmake -S . -B build_2d -DWarpX_DIMS=2 -DWarpX_COMPUTE=OMP \
      -DWarpX_OPENPMD=OFF -DWarpX_QED=OFF -DWarpX_FFT=OFF \
      -DWarpX_PYTHON=OFF -DWarpX_PRECISION=DOUBLE
cmake --build build_2d -j 24

# 1) 上游 ECT 回归（确认 fork 与上游一致）
cd Examples/Tests/embedded_boundary_rotated_cube
OMP_NUM_THREADS=8 ../../../build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB \
    inputs_test_2d_embedded_boundary_rotated_cube
python3 ../../../benchmarks/plotfile_checksum.py diags/diag1000068
#   → 与 Regression/Checksum/benchmarks_json/test_2d_embedded_boundary_rotated_cube.json 对照

# 2) TE（失败）与 TM（成功）对照
cd benchmarks/cases/c02_0_pec
OMP_NUM_THREADS=8 ../../../build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB inputs_2d      # TE: 零/异常
OMP_NUM_THREADS=8 ../../../build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB inputs_2d_tm_flat  # TM: r=-0.9993
python3 dump_tm.py
```

---

## 7. 附录：数据与代码定位

### 7.1 文件清单（本算例目录）

`inputs_2d`（TE，复现失败）、`inputs_2d_tm`（TM，α=30°）、`inputs_2d_tm_flat`（TM，α=π/2，r=−1）、
`inputs_noeb`/`inputs_yee`/`inputs_yee2`（隔离实验）、`analyze.py`/`dump.py`/`dump_tm.py`（分析）；
工具：`../../plotfile_checksum.py`（plotfile checksum）、`../../check_probe.py`（探针检查）、
`../../env_check.py`（环境核对）。

### 7.2 关键代码位置

| 主题 | 位置 |
|---|---|
| XZ 下只算 Rhoy | `EvolveECTRho.cpp`（`EvolveRhoCartesianECT`） |
| ECT B 更新三分支 + Venl | `EvolveB.cpp`（`EvolveBCartesianECT`） |
| E 更新（与 Yee 同码，仅掩膜） | `EvolveE.cpp` |
| 维度限制 abort | `EvolveB.cpp` / `EvolveECTRho.cpp` 的 `WARPX_DIM_3D \|\| WARPX_DIM_XZ` |
| 度量仅 maxLevel | `WarpXInitData.cpp` |

### 7.3 定量结果汇总

rotated_cube：Ex/By 与基准一致（~1e-9）；TE：59 ns 未过峰；TM（α=π/2）：入射 10.02 ns / Ex=1.0000e6，
反射 43.49 ns / Ex=−9.9933e5，**r = −0.9993**（0.07%）。

---

## 8. 结论

1. **不是 fork 回归**：源码逐字节一致 + 上游回归 checksum 对齐（§3.1/§3.2）。
2. **是上游 2D ECT 的实现边界**：XZ 下只实现 TM（仅 `Rhoy`），TE（`Rhox`/`Rhoz`）未实现（§4）。
3. **修正后的基线**：TM + 正入射，PEC 反射 **r = −0.9993**（§3.4c），可作为 C2-0 的合格基线。
4. **后续**：TE 类 2D 算例需转 3D，或补 2D 的 `Rhox`/`Rhoz`（R-B）；完整共形 E 更新
   与对偶几何按评审 D2 单独排期（R-C）。
