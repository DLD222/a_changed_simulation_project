# 宏观介质求解器：输入参考

本文是 `algo.em_solver_medium = macroscopic` 的**用户输入清单**。所有键都在 `macroscopic.*` 下，
除非注明。启动时会打印一段 `[MACRO] ---- effective configuration ----` 摘要，先看它再看结果。

## 1. 最小可用输入

```
algo.em_solver_medium = macroscopic
algo.macroscopic_solver = exponential     # 矩阵指数 M1/M2（唯一支持张量的路径）
algo.macroscopic_subpixel = 1             # 界面子网格平均（有界面时必需）

macroscopic.epsilon  = eps0               # 材料 1
macroscopic.epsilon2 = 4.0*eps0           # 材料 2（有界面时必需）
macroscopic.interface_function(x,y,z) = "-z"   # phi>0 为材料 1，z<0 是材料 1
macroscopic.sigma    = 0.01               # 传导率（标量，全域）
macroscopic.mu       = mu0                # 磁导率（标量，全域）
```

## 2. 全部键

| 键 | 类型 | 说明 |
|---|---|---|
| `epsilon` | 实数 | 材料 1 的介电常数（各向同性） |
| `epsilon_function(x,y,z)` | parser | 空间变化的 ε（材料 1） |
| `epsilon_tensor_xx/yy/zz/xy/xz/yz` | parser ×6 | 材料 1 的 3×3 ε 张量（与上两者**互斥**，同时给会 abort） |
| `epsilon2` | 实数 | 材料 2 的介电常数（有界面时**必需**） |
| `epsilon2_tensor_xx/.../yz` | parser ×6 | 材料 2 的 3×3 ε 张量（与 `epsilon2` 互斥） |
| `sigma` | 实数 | 传导率（各向同性，全域） |
| `sigma_function(x,y,z)` | parser | 空间变化的 σ |
| `sigma_tensor_xx/.../yz` | parser ×6 | 3×3 σ 张量（与 `sigma`/`sigma_function` **互斥**） |
| `mu` | 实数 | 磁导率（各向同性，全域） |
| `mu_function(x,y,z)` | parser | 空间变化的 μ |
| `interface_function(x,y,z)` | parser | 两种材料的分界面，`phi>0` = 材料 1 |

算法开关（在 `algo.*` 下）：

| 键 | 说明 |
|---|---|
| `em_solver_medium = macroscopic` | 打开介质求解器 |
| `macroscopic_solver = exponential` | 矩阵指数 `E^{n+1}=M1 E^n + M2 curlH`，张量必需 |
| `macroscopic_subpixel = 1` | 界面 node 用 τ 平均（各向异性版本） |

## 3. 会被明确拒绝的输入（不再静默忽略）

以下键**会 abort 并给出原因**，避免"设了但没生效"：

```
macroscopic.mu_tensor_*        # μ 张量未实现（内部状态是 B，H=B/μ 只有标量 μ）
macroscopic.sigma2             # 只有一个 σ 场
macroscopic.sigma_function2
```

同时给互相冲突的键也会 abort：

```
macroscopic.epsilon  + epsilon_tensor_*      -> "give either ... not both"
macroscopic.sigma    + sigma_tensor_*        -> 同上
interface_function   但缺 epsilon2           -> "requires a material-2 permittivity"
```

## 4. 已知限制（不是 bug，是未实现）

| 项 | 状态 |
|---|---|
| 材料数 | 只有两种（一个 `interface_function`），不支持三材料/多层 |
| AMR | 介质路径要求 `lev == 0`（多级会 abort） |
| PML | PML 与 `macroscopic` 不兼容（P8） |
| μ 张量 | 未实现 |
| `rho_medium` 诊断 | 该场自动分配、参与 Gauss 检查，但**不是可绘图字段** |
| `sigma_normal_arith` | **内部诊断开关**（默认 0 = 物理解），不要在生产算例里打开 |

## 5. 自检输出

启动时：

```
[MACRO] ---- effective configuration ----
[MACRO]  eps  : 3x3 tensor (material 1)
[MACRO]  eps2 : scalar (material 2)
[MACRO]  sigma: 3x3 tensor
[MACRO]  mu   : isotropic only (tensor mu not implemented)
[MACRO]  interface: -z
[MACRO]  conduction: matrix exponential M1/M2 (tensor capable)
```

运行中（每 `macroscopic.gauss_check_interval` 步，默认见源码）：

```
[GAUSS] step N  max|div D - rho_medium| = ...  max|divD| = ...  max|rho| = ...
```

`max|divD|` 与 `max|rho|` 都为 0 时检查是**空洞**（例如法向入射的 s-pol），不能当作通过。

## 6. 验收基线（法向入射，E3 通量法，带 `R+T` 闭合门）

| 网格 | 基准 `R` | 界面位移 `dz/2` | 良导体 `σ=1e3` |
|---|---|---|---|
| 1024 | 0.1078（−3.0%） | 0.1171（+5.4%） | — |
| 2048 | **0.1113（+0.15%）** | **0.1136（+2.23%）** | **R=1.0000, T=0** |

理论 `R = 0.1111`。**判据必须同时看 `R+T`**（应 ≈1）：否则测量无效。
