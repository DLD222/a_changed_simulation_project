# 多流体矩模型求解器 —— 最终架构（v2 定稿）

<!-- SPDX-License-Identifier: BSD-3-Clause -->


日期：2026-10-31
项目：基于 WarpX 分支（`~/WarpX_fluid`），纯 C++，流体专用（保留粒子能力备用）
目标：n 组分 5/10 矩流体 + **完整麦克斯韦** + **cut cell（必须）** + **电介质 σ/ε/μ** + **AMR 自适应网格（必须）**

---

## 1. 决策汇总（16 项，已逐项确认）

| # | 决策点 | 定稿 |
|---|---|---|
| 1 | 项目形态 | WarpX 全量分支 `~/WarpX_fluid`（不动 `~/WarpX_changed`，源码已复制 2321 文件 / git 已 init） |
| 2 | 场方程 | **完整 Maxwell**：复用 WarpX Yee FDTD + Poisson 散度清洁；GLM 8 变量预留 |
| 3 | 电介质 | σ/ε/μ 三参数；**渐进**：一期常数 → 标量场 → 张量场 → 色散张量场；σ 用 CN 隐式 `1/(1+σdt/2)` |
| 4 | 电介质界面 | **真 cut-cell 嵌入 Riemann**（切向 E/H、法向 D/B 连续，Gkeyll `COPY_B` 思路，自写） |
| 5 | 流体离散 | FV 波传播（Roe/HLLC + MC 限制器），equation-object 预留 DG |
| 6 | 流体模型 | 5 矩 / 10 矩，**非相对论先行**，接口预留 γ 因子（二期相对论） |
| 7 | 十矩闭合 | **k0 正则化（一期）→ 梯度闭合（二期）**，可插拔 |
| 8 | 组分 | **任意 n 组分**，每组分独立 m/q/γ/模型，通用源项与沉积 |
| 9 | 源耦合 | Strang 状态机；一期显式、二期 CN 闭式隐式（ωp/ωc 组合） |
| 10 | cut cell | 必须；**小胞元稳定化 = 三方案原型对比后定稿**（ECT 增强 / 通量重分配 / 局部亚循环） |
| 11 | AMR | 必须；**subcycling**（细网格小 dt + Strang 逐层嵌套）+ reflux + EB 感知保守插值 |
| 12 | 维度/坐标 | **3D + 2D(XZ)** |
| 13 | 精度/GPU | **CPU 双精度优先**，GPU 后期 |
| 14 | EB 几何 | 一期解析隐式函数（parser/level-set），二期 STL 导入 |
| 15 | 粒子 | **保留完整粒子能力**（未来流体+动力学粒子混合） |
| 16 | MVP | **无 AMR 最小闭环**：5 矩单组分 + 完整 Maxwell + cut cell + 常数 σ/ε/μ |

---

## 2. 数学模型（最终）

### 2.1 流体组分 s = 1..n（非相对论）
- **五矩** `U=[ρ, ρux, ρuy, ρuz, E]`，p=(γ−1)(E−½ρ|u|²)，波速 |u|±√(γp/ρ)
- **十矩** `U=[ρ, ρu(3), Pxx_c,Pxy_c,Pxz_c,Pyy_c,Pyz_c,Pzz_c]`，波速 |u|±√(3P11/ρ)、|u|±√(P11/ρ)、u
- 源项：动量 `ρq(E+u×B)`；五矩能量 `J·E`；十矩压强张量磁化旋转 `(q/m)(Ω·P−P·Ω)`
- 闭合（十矩）：一期 k0 超双曲正则化，二期顶点对称化梯度闭合

### 2.2 场（完整 Maxwell + 电介质）
- `∂D/∂t = ∇×H − J`，`∂B/∂t = −∇×E`；`D=εE`、`B=μH`、`J_cond=σE`（σ CN 隐式）
- 散度清洁：`do_dive_cleaning / do_divb_cleaning`（Poisson 投影）
- 介质区无非中性源（J 仅 σE）；界面条件由嵌入 Riemann 保证

### 2.3 耦合（Strang）
```
源(dt/2) → Maxwell 场(dt) → 流体输运(dt) → 源(dt/2)
```
AMR subcycling 下逐层嵌套，跨层 reflux 同步。

---

## 3. 模块架构（`Source/MomentsSolver/`）

```
Source/MomentsSolver/
├── Equations/               Euler5Eqn / TenMomentEqn / MaxwellEqn(GLM预留)  ← 已建接口+Euler5
├── MomentFluidContainer     n 组分 U MultiFab + FV 波传播 + MC 限制器 + 电流/电荷沉积
├── MomentFieldSolver        完整 Maxwell（复用 Yee）+ ε/μ/σ + 散度清洁
├── MomentEMCoupling         Strang 状态机（一期显式 / 二期 CN 隐式）
├── MomentEB                 cut-cell 通量 + 小胞元稳定化(可插拔) + 介质界面嵌入 Riemann
├── MaterialProperties       ε(x)/μ(x)/σ(x)（一期常数，接口预留标量场/张量场/色散）
└── AMR 集成                 subcycling + reflux（复用 WarpX 现有 EB 感知 AMR）
```

复用件：`FiniteDifferenceSolver`（Yee）、`MacroscopicProperties`（材料场组织）、`EmbeddedBoundary/`（EB2、face_areas/edge_lengths、ECT）、`Fluids/MusclHancockUtils.H`（FV 循环骨架）、`HybridPICModel`（输入参数组织）。

---

## 4. 分期路线

| 里程碑 | 内容 | 验收 |
|---|---|---|
| F0 ✅ | 构建 + EB 冒烟 | `warpx.3d.MPI.OMP.DP.PDP.EB` 运行 EB 静电球收敛 |
| F1 ✅ | 方程对象（Euler5/TenMoment/Maxwell）+ FV 波传播驱动 + 单元测试 | 通过：Sod 激波（质量/能量守恒 4e-15，密度轮廓 ρ(0.2/0.5/0.75/0.9)=1.0/0.421/0.265/0.125 与解析一致）；十矩平衡态（100 步零漂移）；真空 EM 传输（脉冲以 c 传播，位移精确 200 格） |
| F2 ✅ | 5 矩单组分 + 完整 Maxwell 显式耦合（无 EB） | 通过：Langmuir 振荡（ωpe 实测 0.99995，误差 5e-5）；两流体 EM 色散（ω_exact=1.4177，实测 1.4171，误差 4.6e-4）。显式源用辛欧拉（对谐振子稳定）+ q/m + 内能守恒能量更新 |
| F3 ✅ | 10 矩 + k0 闭合 | 通过：k0 松弛精确 exp(-νt)、迹守恒；十矩声波 ω 误差 5e-4 |
| M2a（项目 M2） | PEC cut-cell：E+B 共形 + 三方案稳定化对比 → 定稿 | 对比报告 + 导体反射/绕射（已建立 TM 基线 r=−0.9993）<br>**进度（2026-10-31）**：B 侧（2D TE + per-unit-y 几何）完成，L1/L2a/L2b 通过；C 侧（E 场对偶共形 Ampere）与 L3 待做，详见 `ECT_2D_TE_REPORT.md` §14 |
| M2b（项目 M2） | 电介质 σ/ε/μ + 界面嵌入 Riemann | EM 波介质平板 R/T vs Fresnel 解析 |
| M2c（项目 M2） | n 组分任意参数化（m/q/γ/模型） | 多组分耦合回归 |
| M4（项目 M4） | AMR：subcycling + reflux + EB 感知插值 | 自适应算例守恒/精度 |
| 二期 | CN 隐式源、标量/张量场、色散(ADE)、STL、GPU、相对论、DG、GLM | 相应验收 |

> **编号说明（2026-10-31 修正）**：F0–F3 是矩（流体）求解器里程碑；**M2 / M4 沿用项目参考文档**
> （`benchmarks/EM_EB_MHD_VALIDATION_10CASES_v2.md`、`M4_AMR_IMPLEMENTATION_RESEARCH_v1.md`）
> 的编号语义——**M2 = 介质 + EB 单层，M4 = AMR**。此前本文把 cut-cell 叫 M4、介质叫 M5、AMR 叫 M7，
> 与项目编号冲突；现已统一：cut-cell（PEC）= M2a，电介质 = M2b，AMR = M4。

---

## 5. cut-cell 稳定化原型对比（M2a 定稿任务）

| 方案 | 做法 | 指标 |
|---|---|---|
| A ECT 增强 | 复用 `EvolveBCartesianECT`（borrowing/扩展体积）扩展到 E 场+流体 | 局部一阶/全局二阶、侵入性 |
| B 通量重分配/混合 | 小胞元超限通量按体积比分给邻居 | 守恒性、鲁棒性 |
| C 局部亚循环 | 小胞元局部小 dt（AMReX subcycling） | 精度、通信复杂度 |

对比算例：① 平面 EM 波掠斜 PEC（反射精度）；② 超声速流体绕圆柱（守恒误差）；③ 小胞元体积分数扫描（稳定性下限）。

---

## 6. 验证算例套件（建议，待你确认）

- Sod 激波（Euler5 正确性）
- 十矩平衡态 + 各向异性压强波
- Langmuir 振荡（5 矩+Maxwell，ωpe 对标）
- 两流体 EM 波色散（解析对标）
- EM 波垂直入射介质平板：R/T 系数 vs Fresnel 解析（M2b 核心验收）
- 导体球/圆柱绕射（cut cell）
- 小胞元体积分数扫描（稳定性）
- Gkeyll regression 对标（two-fluid shock 等）

---

## 7. 风险与对策

1. **E 场 cut-cell + 介质界面 Riemann（Yee 错位网格）**：全新实现，先做 2D 原型；Meep 法向调和/切向算术公式校核。
2. **AMR subcycling + Strang + cut cell 守恒**：reflux 与体积分数加权插值耦合，M4（AMR）单独攻坚。
3. **十矩 Roe 特征结构**：移植 Gkeyll Maxima 生成代码（`gkyl_wv_ten_moment_priv.h:336`），先单测波分解。
4. **显式源 ωpe 限制**（一期已知）：高密度/磁化算例延后；Strang 接口已预留 CN 隐式。
5. **STL EB 鲁棒性**（二期）：AMReX EB2 STL 已知不稳，需额外几何修复，风险单独管理。

---

## 8. 待确认（开放项）

1. §6 验证算例套件是否满足你的验收预期？是否需补充具体实验/论文对标？
2. 是否按此定稿架构立即开工 F1（方程对象 + FV 驱动 + 单元测试）？

> 历史记录见 `design-plan-moment-solver.md`（含调研笔记与 M0/M1 早期状态）。
