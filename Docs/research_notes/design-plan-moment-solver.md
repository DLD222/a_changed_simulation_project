# 多流体矩模型求解器设计方案（Phase 1）

日期：2026-10-31 · 基于 WarpX（WSL `~/WarpX_changed`）· 调研笔记存档于 artifact 目录 `research-notes-moment-solver.md`

## 0. 用户决策记录

| 决策点 | 选择 |
|---|---|
| 空间离散 | **FV 波传播先行 + 预留 DG 接口**（equation-object 抽象，二期接 modal DG） |
| Cut-cell 小胞元稳定化 | **设计阶段原型对比后定**（ECT 增强 vs 通量重分配/混合 vs 局部亚循环） |
| 电介质能力 | **一期无损 ε/μ + 材料场 MultiFab 接口**（预留空间变化 ε 与二期 σ/Drude-Lorentz 色散） |
| EM 源耦合 | **一期显式 Strang 状态机**（预留二期逐单元 CN 闭式隐式源，Gkeyll `sources_implicit.c` 同款） |

## 1. 数学模型

### 1.1 流体组分 s = 1..N（用户可选模型）
- **五矩**：`U = [ρ, ρu_x, ρu_y, ρu_z, E]`，p = (γ−1)(E − ½ρ|u|²)，通量标准 Euler；波速 |u|±√(γp/ρ)
- **十矩**：`U = [ρ, ρu(3), Pxx_c,Pxy_c,Pxz_c,Pyy_c,Pyz_c,Pzz_c]`；波速 |u|±√(3P11/ρ)（纵）、|u|±√(P11/ρ)（横）、u（4 重熵/剪切波）
- 源项：`ρq(E+u×B)`（动量）、`J·E`（能量，五矩）
- 十矩闭合（按序实现）：① k0 阻尼正则化（Hakim 2008）→ ② 顶点梯度闭合（对称化梯度 + min/max 限制器，注意抛物型 CFL）→ ③ NN 闭合（远期可选）

### 1.2 场（完整麦克斯韦 + 双曲散度清洁）
- 8 变量 `q = [E(3), B(3), φ, ψ]`（Munz/GLM）：φ 传播 div E、ψ 传播 div B，速度 c·e_fact / c·b_fact
- 电介质：ε(x)、μ(x) MultiFab（一期分段常数，介质区无流体源项）；σ 走 CN 隐式 1/(1+σdt/2)（二期）
- 界面条件：切向 E/H 连续（嵌入界面 Riemann 问题，Gkeyll `COPY_B` 思路）；Meep 法向调和/切向算术平均公式作为半胞元等效参数理论参考

### 1.3 耦合（Strang 状态机）
```
源(dt/2) → Maxwell 场(dt) → 流体输运(dt) → 源(dt/2)   [失败回退减 dt 重做，Gkeyll 式]
```
一期源项显式（SSP-RK/前向 Euler）；二期换逐单元 CN 闭式隐式解（ωp·dt、ωc·dt 组合，沿 B 解 3×3）。

## 2. 软件架构（新顶层模块 `Source/MomentsSolver/`）

```
Source/MomentsSolver/
├── Equations/                    ★ equation-object 抽象（DG 二期接口）
│   ├── MomentEquation.H          # 纯虚: num Waves/flux/wave_speeds/Roe 波分解/Riemann 通量
│   ├── Euler5Eqn.{H,cpp}         # 五矩（Roe/HLLC/HLL/Lax 可选）
│   ├── TenMomentEqn.{H,cpp}      # 十矩（特征结构由 Maxima 符号生成，移植 Gkeyll）
│   └── Maxwell8Eqn.{H,cpp}       # GLM Maxwell（Roe 6 波）
├── MomentFluidContainer.{H,cpp}  # 每组分 U MultiFab + 波传播 FV 更新 + MC 限制器
│                                 # （借鉴 Fluids/MusclHancockUtils.H 的循环骨架）
├── MomentFieldSolver.{H,cpp}     # Maxwell 推进（复用 FiniteDifferenceSolver 骨架 + ε/μ 场）
├── MomentStrangDriver.{H,cpp}    # 状态机 + 显式源；二期插 CN 隐式源实现
├── MomentEB.{H,cpp}              # cut-cell：半体积/面积分数通量 + 界面嵌入 Riemann（COPY_B/PEC）
│                                 # + 小胞元稳定化策略可插拔（见 §3 原型对比）
├── MaterialProperties.{H,cpp}    # ε/μ MultiFab（parser/EB 几何驱动；二期 σ/色散）
└── 注册: WarpX.cpp (algo.fluid_solver = moments)、输入参数、Python 绑定、Docs
```

关键复用清单：
- `FiniteDifferenceSolver` Yee 更新骨架（`EvolveE/EvolveB.cpp`）
- `MacroscopicProperties` 的 ε/μ MultiFab 组织方式
- `EmbeddedBoundary/`（EB2 几何、face_areas、edge_lengths、`eb_update_*` 掩膜）
- `Fluids/MusclHancockUtils.H` 的 FV 循环与边界通信骨架
- `HybridPICModel` 的输入参数组织与 Python 暴露方式

## 3. Cut-cell 稳定化原型对比（M4 里程碑）

| 方案 | 原型内容 | 评价指标 |
|---|---|---|
| A. ECT 增强 | 复用 `EvolveBCartesianECT`（borrowing/扩展体积），扩展到 E 与流体通量 | 局部一阶/全局二阶；实现侵入性 |
| B. 通量重分配/混合 | 小胞元超限通量按体积比分给邻居 | 守恒性、鲁棒性、精度损失 |
| C. 局部亚循环 | 小胞元局部小 dt（AMReX subcycling） | 无精度损失；通信复杂度 |

对比算例：①平面 EM 波掠过斜 PEC 表面（反射精度）；②超声速流体绕圆柱（小胞元质量/能量守恒误差）；③ 小胞元体积分数扫描（稳定性下限）。结论写入设计文档后定稿单一方案。

## 4. 里程碑

| 里程碑 | 内容 | 验收 |
|---|---|---|
| M0 | WSL WarpX 构建跑通（含 AMReX_USE_EB） | 现有 regression 一个 EB 算例通过 |
| M1 | Equations/ 三个 equation objects + 单测 | Sod 激波（Euler5）、十矩平衡态微扰、Maxwell8 真空传输 |
| M2 | 5 矩单/双组分 + Maxwell 显式耦合（无 EB） | Langmuir 振荡频率；两流体 EM 色散对标解析 |
| M3 | 10 矩 + k0 闭合（→ 梯度闭合） | 各向异性压强波；Gkeyll ten-moment 回归对标 |
| M4 | cut-cell 三方案原型对比 + 定稿 | §3 对比报告 |
| M5 | 电介质（无损 ε/μ + 界面条件） | EM 波介质平板 R/T vs Fresnel 解析；介质内场解正确 |
| M6 | （二期）CN 隐式源、σ/色散、DG 接口 | 磁化刚性算例；有损介质色散关系 |

## 5. 风险与对策

- **十矩 Roe 特征结构复杂** → 直接移植 Gkeyll Maxima 生成代码（`gkyl_wv_ten_moment_priv.h:336` 起），先做波分解单元测试
- **Yee 错位网格上电介质界面条件** → 界面 Riemann 在面局部系（rot_to_local）求解，参考 Gkeyll `wave_embed_copy_B`；Meep 公式校核
- **显式源 ωpe 限制**（一期已知限制）→ 磁化/高密度算例延后到 M6；Strang 状态机接口已预留
- **梯度闭合抛物型 dt ~ dx²** → 闭合步独立 dt 控制 + 失败回退（Gkeyll 同款 cflrate 机制）

---

## 6. 最终架构决策（已确认 2026-10-31）

| 决策点 | 结论 |
|---|---|
| 语言 | 纯 C++（WarpX 的 C++17/20 + AMReX） |
| 源码策略 | **不动 `~/WarpX_changed`**；整体复制到新目录 `~/WarpX_fluid`（排除 `.git`，已 `git init`） |
| 项目形态 | **A 全量分支**：新增 `Source/MomentsSolver/`，主循环加 `algo.fluid_solver = moments` 分支 |
| 场方程 | **路线 1：复用 WarpX Yee FDTD + Poisson 散度清洁**（`do_dive_cleaning/do_divb_cleaning`）；`MaxwellEqn` 预留 GLM 二期 |
| 粒子处理 | **ii CMake 开关裁剪**：新增 `WarpX_USE_PARTICLES=OFF`，排除动力学粒子/碰撞/QED/Laser 编译 |
| 空间离散 | FV 波传播先行 + equation-object 接口预留 DG |
| Cut-cell 稳定化 | 设计阶段三方案原型对比（ECT 增强 / 通量重分配 / 亚循环）后定稿 |
| 电介质 | 一期无损 ε/μ + 材料场 MultiFab；二期 σ/Drude-Lorentz |
| 源耦合 | 一期显式 Strang 状态机；二期 CN 闭式隐式源 |

**项目路径**：`~/WarpX_fluid`（`//wsl.localhost/Ubuntu-test<repo>`）
**基线构建配置**：`cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=OMP -DWarpX_OPENPMD=OFF -DWarpX_QED=OFF -DWarpX_FFT=OFF -DWarpX_PYTHON=OFF -DWarpX_PRECISION=DOUBLE`（EB 默认 ON，MPI ON）

**M0 状态**：✅ 完成。源码复制（2321 文件 / 23M / git 已 init）→ 构建成功（`warpx.3d.MPI.OMP.DP.PDP.EB`，3D + EB + 双精度 + OMP/MPI）→ EB 冒烟测试通过（3D 静电球 EB，MLMG cut-cell Poisson 收敛，正常 `finalized`）。

**M1 状态**：🚧 已建 `Source/MomentsSolver/Equations/` 骨架（`MomentEquation.H` 方程接口 + `Euler5Eqn.H` 五矩方程对象，vtable 函数指针模式，GPU 友好）。待接 FV 驱动与单元测试。
