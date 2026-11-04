# 场 / 等离子体 / 介质 模块化架构建议（建议稿，非决定）

> 本文是**建议稿**：给出目标、可选实现路径、扩展配方与注意事项，供决策参考。
> 文中标注“建议 / 可选 / 待定”的地方**不是拍板**；最终取舍由项目负责人决定。
> 日期：2026-11-04 · 基于 `~/WarpX_fluid` 源码阅读与公开文献。

---

## 0. 文档定位

- **目的**：为“场、等离子体、介质分块 + 极致可扩展”提供一份可讨论的架构方案。
- **不是**：实现承诺、时间表、或对现有代码的否定。
- **读者**：项目负责人、后续实现者。
- **原则**：不推倒 AMReX/WarpX；以适配器与 strangler 方式渐进；每个阶段保持现有回归。

---

## 1. 目的与范围

### 1.1 要解决的问题

当前扩展一个新物理（例如“粒子间引力”“介质相变/气化”）通常需要同时改动：

```
WarpX.H/.cpp（加成员、加调用）
Fields.H（加 FieldType 枚举）
Evolve/WarpXEvolve.cpp（加时间步顺序）
FieldSolver/（新求解器）
Particles/ 或 Fluids/（新沉积/推进）
Initialization/（新输入解析）
Diagnostics/（新输出）
```

即“一处物理、多处改动”，违反开闭原则，也让回归面变大。

### 1.2 目标（建议）

1. **统一等离子体块**：粒子与流体（含 5/10 矩）在同一模块内，共享物种、EOS、闭合、摩擦/碰撞、沉积/收集接口。
2. **场块插件化**：Maxwell（Yee/PSATD/CKC）、Poisson（静电/引力）、其他场方程都是可替换插件。
3. **介质块独立**：ε/σ/μ、子网格/界面、面电荷、色散（ADE）、**相变**（气化/沸腾）都是可替换的材料模型。
4. **耦合层显式**：相互作用（电磁、引力、碰撞/摩擦）是独立单元，有明确的 deposit/solve/gather 端口。
5. **扩展成本**：新增一个“场 + 力”的相互作用，理想情况下只写 **2 个 .cpp**（一个场方程 + 一个相互作用），加注册与输入；新增一个材料过程（相变），理想情况下只写 **2 个 .cpp**（一个相变模型 + 一个物性模型）。
6. **可测试**：每个端口有 mock；新模型带契约测试；现有 regression 逐位/数值不变。

### 1.3 非目标（建议）

- 不重写 AMReX 的网格/并行/GPU 基础设施；
- 不在第一阶段引入 DG/高阶（保留接口即可）；
- 不追求一次到位；允许分阶段、允许旧路径暂时共存。

---

## 2. 现状审计（代码证据）

### 2.1 结构

```
Source/
├─ WarpX.H / WarpX.cpp          # 上帝类：参数、场、粒子、流体、诊断、EB、PML、时间循环
├─ Fields.H                     # FieldType 平铺枚举（~100 项）
├─ FieldSolver/
│  ├─ FiniteDifferenceSolver/   # Yee/CKC + MacroscopicProperties（介质嵌在这里）
│  ├─ SpectralSolver/           # PSATD
│  ├─ ElectrostaticSolvers/     # 静电
│  ├─ ImplicitSolvers/          # 隐式
│  ├─ WarpXPushFieldsEM.cpp     # EvolveE/B + 介质 ρ_medium + Gauss 诊断
│  └─ ...
├─ Particles/                   # PIC
├─ Fluids/                      # 冷相对论流体（N, NU；无压强/摩擦）
├─ MomentsSolver/               # 1D 5/10 矩标准测试床（Euler5、FV、耦合、单元测试）
├─ Evolve/WarpXEvolve.cpp       # 时间循环与耦合顺序
├─ Diagnostics/ · BoundaryConditions/ · Initialization/ · Python/ · ablastr/
```

### 2.2 关键耦合点（读到的）

| 位置 | 作用 | 问题 |
|---|---|---|
| `WarpX` 类 | 持有全部状态与求解器 | 上帝对象；扩展需改类 |
| `FieldType` 枚举 | 场注册表 | 平铺、语义混杂（状态/源/临时/EB/PML） |
| `Evolve()` | 时间循环 + 耦合顺序 | 物理顺序写死在核心 |
| `MacroscopicProperties` | 介质 ε/σ + 子网格 + ρ_medium | 嵌在 FDTD 目录，非独立块 |
| `Fluids/` 与 `MomentsSolver/` | 两套流体 | 不共享 EOS/闭合/摩擦接口 |
| `MultiFabRegister` | 场容器 | 好用，但缺少类型化语义 |

### 2.3 可复用资产（建议保留）

- `ablastr::fields::MultiFabRegister`：场容器，适合做“状态注册表”的底座；
- `FieldSolver` 的 Yee/CKC/PSATD 实现：作为场方程插件的**适配对象**；
- `MomentsSolver/Equations/` 的 equation-object（Euler5/TenMoment/Maxwell6）：已是“物理即插件”的雏形；
- `MomentsSolver/Coupling/` 的 Strang 状态机：可扩展为耦合层参考；
- 现有 regression（P0/E3/S5/E6/L4/L5/L6…）：迁移的安全网。

---

## 3. 文献调研（架构模式）

| 来源 | 机制 | 对本项目的启示 |
|---|---|---|
| **MOOSE** | `Kernel`、`Material`、`InterfaceKernel` 插件；输入卡写类型名；按物理模块注册（[概览](https://icmed.engin.umich.edu/wp-content/uploads/sites/176/2018/05/Summer-School-MOOSE-talk-2018-reduced.pdf)） | 场/材料/界面都是“对象 + 注册 + 输入” |
| **Smilei** | 物种 / 场求解器 / patch 分层；模块化 PIC（[论文](https://arxiv.org/html/1702.05128v1)） | 粒子与场解耦，interaction 独立 |
| **PIConGPU** | species 属性 + 插件式 current deposition（[仓库](https://github.com/computationalradiationphysics/picongpu)） | 沉积/收集是端口，可换 |
| **Gkeyll** | equation-object（flux/wave speed/Riemann）+ 波传播驱动 | 模型即插件；一维测试床先行 |
| **相变/沸腾** | 相场（Cahn–Hilliard/Allen–Cahn）、焓法、混合物介电（Maxwell–Garnett/Bruggeman）、潜热；[电场诱导相分离](https://pubs.aip.org/aip/jcp/article/161/4/044704/3304666/Electric-field-induced-phase-separation-and)、[C6F12O 沸腾](https://ouci.dntb.gov.ua/en/works/42B0kXV4/) | 介质相变 = 相变模型 + 物性模型 |
| **自引力** | 粒子-网格（PM）Poisson：质量沉积 + 力 gather；[数值宇宙学讲义](https://indico.iap.fr/event/25/contributions/223/attachments/153/278/LesHouches_Lecture1.pdf) | 引力 = 场方程 + 相互作用 |

**共同模式**：

```
物理 = 状态 + 模型（插件）+ 端口
耦合 = 显式的 interaction/coupling 层
扩展 = 实现接口 + 注册 + 输入
```

---

## 4. 目标架构（建议）

### 4.1 分层图

```
┌──────────────────────────── Core / Runtime ────────────────────────────┐
│  typed state registry · geometry/AMR · clock/time-loop · diagnostics   │
│  particle/fluid containers · memory/GPU · Python/inputs                │
└───────────────▲────────────────▲────────────────▲─────────────────────┘
                │ ports          │ ports          │ ports
┌───────────────┴──────┐ ┌───────┴────────┐ ┌────┴──────────────────────┐
│   Fields block       │ │  Plasma block  │ │   Medium block            │
│  E,B / φ_g / ...     │ │ species (PIC,  │ │ phase, T, ε, σ, μ         │
│  Maxwell, Poisson    │ │ fluid, moment) │ │ transition model          │
│  boundary/PML        │ │ EOS, closure   │ │ property model            │
│  source ports        │ │ interactions   │ │ interface/subpixel/charge │
└───────────────▲──────┘ └───────▲────────┘ └────▲──────────────────────┘
                └────────── Coupling layer ───────┘
                interaction (EM, gravity, collisions/friction) + split
```

### 4.2 核心接口（建议稿）

```cpp
// 等离子体：粒子/流体/矩统一
class ISpecies {
public:
    virtual void push(const FieldView&, Real dt) = 0;
    virtual void deposit(FieldSource&) = 0;      // J, rho, rho_m
    virtual void source(const FieldView&, Real dt) = 0;
    virtual void transport(Real dt) = 0;         // 流体/矩输运
    virtual ~ISpecies() = default;
};

// 场方程：Maxwell / Poisson / ...
class IFieldEquation {
public:
    virtual void solve(FieldState&, const FieldSource&, Real dt) = 0;
    virtual void boundary(FieldState&) = 0;
    virtual ~IFieldEquation() = default;
};

// 相互作用：电磁 / 引力 / 碰撞
class IInteraction {
public:
    virtual void deposit(const PlasmaState&, FieldSource&) = 0;
    virtual void solve(FieldState&, const FieldSource&, Real dt) = 0;
    virtual void gather(PlasmaState&, const FieldState&, Real dt) = 0;
    virtual ~IInteraction() = default;
};

// EOS / 闭合 / 摩擦
class IEquationOfState {
public:
    virtual Real p(Real rho, Real e) const = 0;
    virtual Real cs2(Real rho, Real e) const = 0;
    virtual Real T(Real rho, Real e) const = 0;
    virtual ~IEquationOfState() = default;
};
class IClosure { public: virtual void apply(SpeciesState&, Real dt) = 0; };
class IFriction { public: virtual void apply(std::span<SpeciesState>, Real dt) = 0; };

// 介质：物性 / 相变 / 响应
class IMediumProperty {
public:
    virtual void fill(Level&, const MediumState&) = 0;
    virtual ~IMediumProperty() = default;
};
class IPhaseTransition {
public:
    virtual void advance(MediumState&, const FieldView&, Real dt) = 0;
    virtual ~IPhaseTransition() = default;
};
class IMediumResponse {
public:
    virtual void respond(const FieldView&, FieldSource&, Real dt) = 0;
    virtual ~IMediumResponse() = default;
};
```

### 4.3 状态注册表（建议）

- 把 `FieldType` 平铺枚举改为**分层的 typed key**：

```cpp
enum class FieldRole { State, Source, Aux, Scratch };
template <FieldRole R, class T> struct FieldKey { ... }; // 编译期语义
```

- `FieldState` 提供按 key 的访问；`phi_g`、`T_medium`、`phase` 等新状态以 open slot 注册，**不改核心枚举**。

### 4.4 注册与输入（建议）

```cpp
register_interaction("gravity", ...);
register_field_equation("gravity", ...);
register_phase_transition("vaporization", ...);
register_medium_property("dielectric_mixture", ...);
```

```
interactions = em gravity
gravity.poisson_solver = fft
medium.model = dielectric_mixture
medium.transition = vaporization
```

---

## 5. 扩展配方（“两个 .cpp”的例子）

### 5.1 例 1：粒子间引力

**新增文件（建议）**

```
Gravity/.cpp:
  GravityFieldEquation.cpp   // 质量密度 → Poisson 解 φ_g
  GravityInteraction.cpp     // deposit ρ_m；gather -∇φ_g；push
```

**不改**：WarpX 核心、Evolve 顺序、FieldType 枚举（`phi_g` 走 open slot）、诊断（新状态自动可输出）。

**注册与输入**：见 §4.4。

**验证**：均匀质量分布的解析势；双粒子自由落体；总能量（动能+引力势）守恒；与 PM 参考对比。

### 5.2 例 2：介质相变（气化/沸腾）

**新增文件（建议）**

```
Medium/PhaseTransition/:
  VaporizationModel.cpp      // T、潜热、相分数、焓法/相场速率
  DielectricMixture.cpp      // ε(phase,T)、σ(phase,T)、界面 subpixel
```

**不改**：场求解器（`IMediumResponse` 把 σE/ρ_medium 喂给源端口）、时间循环。

**物理选项（待定）**：

- 尖锐界面：焓法 + VOF/level-set；
- 扩散界面：Cahn–Hilliard / Allen–Cahn 相场；
- 介电混合：Maxwell–Garnett / Bruggeman / 线性插值；
- 潜热/能量：与流体能量方程或介质能量方程耦合。

**验证**：一维 Stefan 问题；均匀场下相变前沿速度；能量守恒；与解析/文献对比。

### 5.3 例 3（可选）：组分间摩擦

```
Friction/:
  SimpleDrag.cpp      // 常数 ν，精确指数
  Braginskii.cpp      // ν ∝ n Z² lnΛ / T^{3/2}，冻结 ν 指数/CN
  UserDefined.cpp     // 外部回调
```

**不改**：5 矩的输运/场耦合；摩擦作为耦合层的 source。

---

## 6. 实现方法（建议分阶段）

### P0：接口 + mock + 适配器

```
- 定义 §4.2 的接口；
- 为现有 FiniteDifferenceSolver/SpectralSolver/MultiParticleContainer/WarpXFluidContainer/
  MacroscopicProperties 写适配器；
- Evolve() 改为“适配器 + Coupling”调用；
- 验收：现有 regression 逐位/数值不变。
```

### P1：Plasma 块

```
- EOS 接口 + 理想气体/等温/用户；
- 5 矩（先在 MomentsSolver 1D 测试床）→ 适配器；
- 摩擦三模型（简单/Braginskii/用户），指数积分；
- 验证：Sod、声波、双组分摩擦弛豫、动量/能量守恒。
```

### P2：Fields 块

```
- FieldState/FieldSource + typed registry；
- 场方程插件（Yee/PSATD/CKC/Poisson）；
- 验收：P0/E3/S5/E6 逐位。
```

### P3：Medium 块

```
- 物性/相变/响应接口；
- ε/σ/μ 张量 + subpixel + 界面 + ρ_medium；
- 相变示例；
- 验收：L5/L6 + 相变解析。
```

### P4：Coupling 层

```
- interaction 注册与顺序；
- Strang/IMEX；
- 能量/电荷审计。
```

### P5：strangler 淘汰旧路径

```
- 旧 FieldType/旧 Evolve 分支逐步删除；
- 文档与测试同步。
```

### 数据/性能/测试（建议）

```
- 数据：SoA/AMReX；状态归 Core，块只拿 view；GPU 热路径静态多态；
- 性能：虚函数只在每步一次；kernel 内联；避免 std::function；
- 测试：每个端口 mock；新模型契约测试；1D 先行；2D/3D 回归；
- 兼容：新老 key 双写；Python 输入逐步迁移。
```

---

## 7. 注意事项与风险（建议逐条评估）

### 7.1 性能

- 抽象层可能拖慢 GPU 热路径 → 建议静态多态 + 内联，虚函数只在“每步一次”的层；
- 状态注册表的间接访问 → 建议 view/span，kernel 内只读裸指针。

### 7.2 AMReX/GPU

- MultiFab 生命周期与 view 的边界 → 状态归 Core，块只拿 view；
- 新状态（φ_g、T、phase）的 ghost/边界/AMR 同步 → 需要统一的注册钩子。

### 7.3 迁移

- 现有 Evolve 顺序敏感 → 适配器阶段保留原顺序；
- `FieldType` 平铺枚举迁移成本 → 双写 + 编译期映射；
- 回归覆盖面 → 每个阶段都要 P0/E3/S5/E6/L4/L5/L6 逐位/数值。

### 7.4 物理正确性

- 相变的热力学闭合（潜热、焓、相场）需要专门验证；
- 引力的 PM 自引力、边界、周期 vs 开放、软化；
- 摩擦的刚性、动量/能量守恒、熵产；
- 介质相变后的 ε/σ 与 EM 界面条件（子网格/面电荷）。

### 7.5 过度设计

- 只为**已知扩展点**（EM/引力/碰撞/相变/EOS）建接口；
- 不为假想需求建抽象；每个抽象配一个真实用例。

### 7.6 组织与流程

- 接口变更需要评审；
- 新物理 = 接口 + 注册 + 契约测试 + 文档；
- 迁移期间保持“旧路径可回退”。

---

## 8. 参考

1. MOOSE 框架概览：<https://icmed.engin.umich.edu/wp-content/uploads/sites/176/2018/05/Summer-School-MOOSE-talk-2018-reduced.pdf>
2. Smilei：<https://arxiv.org/html/1702.05128v1>
3. PIConGPU：<https://github.com/computationalradiationphysics/picongpu>
4. Gkeyll 两流体（5/10 矩）：<https://ar5iv.labs.arxiv.org/html/1903.09618>
5. 电场诱导相分离（相场/介电混合）：<https://pubs.aip.org/aip/jcp/article/161/4/044704/3304666/Electric-field-induced-phase-separation-and>
6. C6F12O 沸腾/相变：<https://ouci.dntb.gov.ua/en/works/42B0kXV4/>
7. 数值宇宙学（PM 引力）：<https://indico.iap.fr/event/25/contributions/223/attachments/153/278/LesHouches_Lecture1.pdf>
8. WarpX 概览：<https://indico.slac.stanford.edu/event/8597/contributions/8188/attachments/3935/10601/2024_beam_beam.pdf>

---

## 9. 建议的决策清单（留给负责人，不在本文拍板）

| # | 决策点 | 选项 | 备注 |
|---|---|---|---|
| 1 | 状态注册表 | 保留 `FieldType` / 改 typed registry / 双写过渡 | 影响编译期语义与迁移成本 |
| 2 | 多态方式 | 静态多态优先 / 运行时虚函数 / 混合 | 影响性能与扩展便利 |
| 3 | 统一等离子体块 | 粒子+流体同块 / 分块但共享接口 | 影响物种/EOS/摩擦的复用 |
| 4 | interaction 粒度 | 场+力一体 / 场与力分离 | 引力/EM/碰撞的复用方式 |
| 5 | 相变模型 | 尖锐（焓/VOF）/ 扩散（相场）/ 两者 | 影响界面与 EM 耦合 |
| 6 | 引力求解器 | FFT（周期）/ MLMG（开放）/ 两者 | 影响边界与精度 |
| 7 | 摩擦模型 | 简单 / Braginskii / 用户 | 影响刚性与守恒 |
| 8 | 迁移范围 | 全量 / 只新物理走新接口 | 影响回归风险 |
| 9 | 1D 测试床与 3D 的关系 | 先 1D 后 3D / 直接 3D | 影响验证速度 |
| 10 | 旧路径淘汰时机 | 立即 / 迁移完成后 | 影响维护成本 |

---

## 10. 一句话

**本报告建议：以“Core + Plasma（粒子/流体同块）+ Fields（场方程插件）+ Medium（物性/相变插件）+ Coupling（interaction）”为目标，用 typed 状态、工厂注册、静态多态优先与契约测试支撑扩展；引力与相变各用“场方程+相互作用”“相变模型+物性模型”两个文件即可接入；迁移用适配器 + strangler + 逐位回归；所有取舍（注册表、多态、相变模型、引力求解器、迁移范围等）列在 §9，交由负责人决定，本文不拍板。**
