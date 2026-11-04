# 场 / 等离子体 / 介质 模块化与通用可扩展架构建议（建议稿，非决定）

> 本文是**建议稿**：给出目标、通用扩展机制、可选实现路径、风险与决策点。
> 文中“建议 / 可选 / 待定”均**不是拍板**；最终取舍由项目负责人决定。
> 说明：文中的“引力”“相变”等只是**成本量级**的示意，**不是**架构针对的对象；
> 架构面向的是“任意新场 / 任意新耦合 / 任意新物种模型 / 任意新材料过程”。
> 日期：2026-11-04 · 基于 `~/WarpX_fluid` 源码阅读与公开文献。

---

## 0. 文档定位

- **目的**：为“场、等离子体、介质分块 + 通用可扩展”提供可讨论的方案。
- **不是**：实现承诺、时间表、对现有代码的否定，或针对某两个示例的专门设计。
- **原则**：不推倒 AMReX/WarpX；适配器 + strangler 渐进；每阶段保持现有回归。

---

## 1. 目的：通用可扩展，而不是“给两个例子开口子”

### 1.1 问题

当前扩展任何新物理（不限于引力/相变）通常要同时改：

```
WarpX.H/.cpp（加成员、加调用）
Fields.H（加 FieldType 枚举）
Evolve/WarpXEvolve.cpp（加时间步顺序）
FieldSolver/（新求解器）
Particles/ 或 Fluids/（新沉积/推进）
Initialization/（新输入解析）
Diagnostics/（新输出）
```

即“一处物理、多处改动”，扩展成本随需求种类线性增长，回归面也随之扩大。

### 1.2 目标（建议）

1. **扩展点分类**：把“新物理”归入有限个**通用角色**（见 §2），每个角色有稳定的接口。
2. **任意组合**：新物理可以是“场 + 耦合”“物种模型 + 闭合”“材料过程 + 物性”“约束 + 源”，也可以是它们的组合；架构不预设具体物理。
3. **声明式接入**：插件声明自己**读/写哪些状态、需要哪些端口、以什么顺序、是否守恒**；调度器按声明组合。
4. **成本可控**：新增一个角色通常 1–2 个文件 + 注册 + 输入 + 契约测试；不碰核心。
5. **可测试**：每个端口有 mock；新插件带契约测试；现有 regression 逐位/数值不变。

### 1.3 非目标（建议）

- 不重写 AMReX 的网格/并行/GPU；
- 不在第一阶段引入 DG/高阶（保留接口）；
- 不追求一次到位；允许旧路径暂时共存。

---

## 2. 通用扩展点分类（架构的核心）

| 扩展点 | 通用接口 | 可容纳的“各式各样”需求 | 新增文件量级 |
|---|---|---|---|
| **场方程** | `IFieldEquation` | 任意椭圆/抛物/双曲场：Maxwell、Poisson、输运、扩散、约束场、频域/时域、线性/非线性 | 1 |
| **相互作用/耦合** | `IInteraction` | 任意耦合：局部/全局、成对/场介导、碰撞/摩擦、加热/电离/复合、外力 | 1 |
| **物种模型** | `ISpeciesModel` | 粒子/流体/矩/混合；任意闭合；任意状态变量 | 1 |
| **材料模型** | `IMaterialModel` | 物性、状态演化、相变、界面、损伤、化学反应、多相混合 | 1–2 |
| **约束** | `IConstraint` | Gauss、divB、不可压、准中性、总能量/总动量约束 | 1 |
| **源/汇** | `ISource` | 外源、加热、电离、复合、注入、吸收 | 1 |
| **边界** | `IBoundary` | 新 BC、PML 变体、界面条件 | 1 |
| **诊断** | `IDiagnostic` | 新输出、新分析、新误差估计 | 1 |
| **时间积分/分裂** | `ITimeIntegrator` | 新 splitting、IMEX、子循环、自适应 | 1 |
| **闭合/EOS/摩擦** | `IClosure` / `IEquationOfState` / `IFriction` | 材料/动理学闭合、状态方程、碰撞模型 | 1 |

**关键**：这些是**角色**，不是具体物理。任何新需求先问“它属于哪个角色/哪几个角色的组合”，而不是“给它单独开一个分支”。

---

## 3. 现状审计（代码证据）

### 3.1 结构

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
├─ MomentsSolver/               # 1D 5/10 矩标准测试床（equation-object、FV、耦合、单测）
├─ Evolve/WarpXEvolve.cpp       # 时间循环与耦合顺序
├─ Diagnostics/ · BoundaryConditions/ · Initialization/ · Python/ · ablastr/
```

### 3.2 关键耦合点

| 位置 | 作用 | 问题 |
|---|---|---|
| `WarpX` 类 | 持有全部状态与求解器 | 上帝对象；扩展需改类 |
| `FieldType` 枚举 | 场注册表 | 平铺、语义混杂（状态/源/临时/EB/PML） |
| `Evolve()` | 时间循环 + 耦合顺序 | 物理顺序写死在核心 |
| `MacroscopicProperties` | 介质 ε/σ + 子网格 + ρ_medium | 嵌在 FDTD 目录，非独立块 |
| `Fluids/` 与 `MomentsSolver/` | 两套流体 | 不共享 EOS/闭合/摩擦接口 |
| `MultiFabRegister` | 场容器 | 好用，但缺类型化语义 |

### 3.3 可复用资产（建议保留）

- `ablastr::fields::MultiFabRegister`：状态注册表的底座；
- `FieldSolver` 的 Yee/CKC/PSATD：场方程插件的适配对象；
- `MomentsSolver/Equations/` 的 equation-object：已是“物理即插件”的雏形；
- `MomentsSolver/Coupling/` 的 Strang：耦合层参考；
- 现有 regression（P0/E3/S5/E6/L4/L5/L6…）：迁移安全网。

---

## 4. 文献调研（通用模式）

| 来源 | 机制 | 启示 |
|---|---|---|
| **MOOSE** | `Kernel`、`Material`、`InterfaceKernel` 插件；输入写类型名 | 通用角色 + 注册 + 输入（[概览](https://icmed.engin.umich.edu/wp-content/uploads/sites/176/2018/05/Summer-School-MOOSE-talk-2018-reduced.pdf)） |
| **Smilei** | 物种/场求解器/patch 分层 | 粒子与场解耦，interaction 独立（[论文](https://arxiv.org/html/1702.05128v1)） |
| **PIConGPU** | species 属性 + 插件式沉积 | 沉积/收集是端口（[仓库](https://github.com/computationalradiationphysics/picongpu)） |
| **Gkeyll** | equation-object（flux/wave speed/Riemann） | 模型即插件 |
| **相场/多相** | Cahn–Hilliard/Allen–Cahn、焓法、混合物介电 | 材料过程 = 状态 + 速率 + 物性（[相分离](https://pubs.aip.org/aip/jcp/article/161/4/044704/3304666/Electric-field-induced-phase-separation-and)） |
| **PM 引力** | 质量沉积 + Poisson + 力 gather | 场+耦合的模式（[讲义](https://indico.iap.fr/event/25/contributions/223/attachments/153/278/LesHouches_Lecture1.pdf)） |

**共同模式**：物理 = 状态 + 模型插件 + 端口；耦合 = 显式层；扩展 = 接口 + 注册 + 输入 + 契约测试。

---

## 5. 目标架构（建议）

### 5.1 分层图

```
┌──────────────────────────── Core / Runtime ────────────────────────────┐
│  typed state registry · geometry/AMR · clock · scheduler · diagnostics │
│  containers · memory/GPU · Python/inputs                               │
└───────────────▲────────────────▲────────────────▲─────────────────────┘
                │ ports          │ ports          │ ports
┌───────────────┴──────┐ ┌───────┴────────┐ ┌────┴──────────────────────┐
│   Fields block       │ │  Plasma block  │ │   Medium block            │
│  任意场状态/方程      │ │ 任意物种模型   │ │  任意材料状态/过程        │
│  边界/约束/源端口     │ │ EOS/闭合/耦合  │ │  物性/界面/响应           │
└───────────────▲──────┘ └───────▲────────┘ └────▲──────────────────────┘
                └────────── Coupling layer ───────┘
                interactions + constraints + sources + split/IMEX
```

### 5.2 通用接口（建议稿）

```cpp
// 场：任意场方程
class IFieldEquation {
public:
    virtual void solve(FieldState&, const FieldSource&, Real dt) = 0;
    virtual void boundary(FieldState&) = 0;
    virtual ~IFieldEquation() = default;
};

// 物种：粒子/流体/矩/混合
class ISpeciesModel {
public:
    // advance = PIC 的粒子推 / 流体的输运 / 矩的输运（同一个语义）
    virtual void advance(const FieldView&, Real dt) = 0;
    virtual void deposit(FieldSource&) = 0;      // 任意源：J, rho, rho_m, ...
    virtual ~ISpeciesModel() = default;
};

// 耦合：任意相互作用
class IInteraction {
public:
    // 相互作用只负责“源/力”两侧；场方程由它持有/引用的 IFieldEquation 完成
    virtual void deposit(const PlasmaState&, FieldSource&) = 0;
    virtual void gather(PlasmaState&, const FieldState&, Real dt) = 0;
    // 可选：声明它使用哪个场方程/状态端口
    virtual FieldEquationId field_equation() const = 0;
    virtual ~IInteraction() = default;
};

// 材料：物性 / 状态过程 / 响应（三者消费者不同，建议分开）
class IMaterialProperty {   // 场求解器要的 ε/σ/μ（可缓存）
public:
    virtual void fill(Level&, const MaterialState&) = 0;
    virtual ~IMaterialProperty() = default;
};
class IMaterialProcess {    // 相变/反应/损伤等状态演化
public:
    virtual void advance(MaterialState&, const FieldView&, Real dt) = 0;
    virtual ~IMaterialProcess() = default;
};
class IMaterialResponse {   // 源：σE、ρ_medium、表面电荷
public:
    virtual void respond(const FieldView&, FieldSource&, Real dt) = 0;
    virtual ~IMaterialResponse() = default;
};

// 约束 / 源 / 边界 / 诊断 / 积分器
class IConstraint  { public: virtual void project(State&) = 0; };
class ISource      { public: virtual void apply(State&, Real dt) = 0; };
class IBoundary    { public: virtual void apply(State&) = 0; };
class IDiagnostic  { public: virtual void evaluate(const State&) = 0; };
class ITimeIntegrator { public: virtual void step(State&, Real dt) = 0; };
```

### 5.3 插件的“声明”与调度

```cpp
struct PluginSpec {
    std::vector<StateKey> reads;      // 读哪些状态
    std::vector<StateKey> writes;     // 写哪些状态
    std::vector<Port>     ports;      // 需要哪些端口
    int                   order;      // 调度顺序
    bool                  conservative; // 是否需守恒审计
    std::optional<Real>   dt_scale;   // 子循环/分裂
};
```

调度器按 `reads/writes/order` 组成**同一 split 阶段内的偏序**；**场与源本质上是循环依赖**（场依赖源、源依赖场），
不能靠 reads/writes 的 DAG 排序，必须由 operator split（Strang）、IMEX 或定点迭代显式处理循环。
插件之间**只通过状态端口**交互，不互相依赖。

### 5.4 状态注册表（建议）

- 把 `FieldType` 平铺枚举改为 **typed key**：

```cpp
enum class FieldRole { State, Source, Aux, Scratch };
template <FieldRole R, class T> struct FieldKey { ... };
```

- 新状态（任意名称/类型）以 open slot 注册，**不改核心枚举**；`FieldState` 按 key 访问。
- 注意：编译期 `FieldKey<R,T>` 只能覆盖**已知**状态；插件自定义的**运行时命名 slot**
  （string/id + 布局）需要类型擦除的 registry，二者要分清。

### 5.5 注册与输入（建议）

```cpp
register_field_equation("...", ...);
register_interaction("...", ...);
register_species_model("...", ...);
register_material_model("...", ...);
register_constraint("...", ...);
```

```
fields.equations = maxwell ...
interactions = em ...
species.models = pic ...
material.models = ...
```

输入是**通用的**：插件自己解析自己的段；核心只负责把名字映射到工厂。

---

## 6. 通用扩展配方（“新增一个东西”的步骤）

> 以下是**角色**层面的步骤；具体物理不限。

```
1. 判定角色：场方程 / 相互作用 / 物种模型 / 材料模型 / 约束 / 源 / 边界 / 诊断 / 积分器
2. 实现接口（1 个 .cpp，复杂材料可拆 2 个：状态过程 + 物性）
3. 声明 PluginSpec（reads/writes/ports/order/dt_scale/conservative）
4. 注册到对应工厂
5. 加输入段（插件自己的参数）
6. 加状态 slot（若需要新状态）
7. 契约测试（mock + 1D）+ 一个真实算例
8. 文档
```

**成本量级**：单角色 **1–2 个核心 .cpp** + 状态 slot/边界/输入/诊断/契约测试等“胶水”；
“不碰核心”指不碰 WarpX/Evolve/FieldType 的分支逻辑，但**胶水不可省**。

**组合示例（仅示组合能力，不是设计对象）**：

```
- 一个“场方程 + 相互作用” = 新力场（任意：引力/静磁/约束场/…）
- 一个“材料状态过程 + 物性” = 新材料过程（任意：相变/损伤/化学反应/…）
- 一个“物种模型 + 闭合” = 新流体/矩模型
- 一个“约束 + 源” = 新守恒/清洁方案
- 以上可叠加：多场、多相互作用、多材料过程、链式过程
```

**关键**：这些组合由 **PluginSpec/调度器** 自动编排，而不是为每个组合写死分支。

---

## 7. 实现方法（建议分阶段）

### P0：接口 + mock + 适配器

```
- 定义 §5.2 接口与 §5.3 PluginSpec；
- 为现有 FiniteDifferenceSolver/SpectralSolver/MultiParticleContainer/WarpXFluidContainer/
  MacroscopicProperties 写适配器；
- Evolve() 改为“适配器 + 调度器”；
- 验收：现有 regression 逐位/数值不变。
```

### P1：Plasma 块

```
- EOS/闭合/摩擦接口（先在 MomentsSolver 1D 测试床验证）；
- 5 矩/10 矩/冷流体统一到 ISpeciesModel；
- 验证：Sod、声波、双组分弛豫、守恒。
```

### P2：Fields 块

```
- FieldState/FieldSource + typed registry；
- 场方程插件（Yee/PSATD/CKC/Poisson/…）；
- 验收：P0/E3/S5/E6 逐位。
```

### P3：Medium 块

```
- 材料状态/过程/物性/界面/响应接口；
- ε/σ/μ 张量 + subpixel + 面电荷；
- 材料过程示例；
- 验收：L5/L6 + 解析。
```

### P4：Coupling 层

```
- interaction/constraint/source 注册与 DAG 调度；
- Strang/IMEX；能量/电荷/动量审计。
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
- 测试：每个端口 mock；新插件契约测试；1D 先行；2D/3D 回归；
- 兼容：新老 key 双写；Python 输入逐步迁移。
```

---

## 8. 注意事项与风险（逐条评估）

### 8.1 性能

- 抽象层可能拖慢 GPU 热路径 → 静态多态 + 内联，虚函数只在“每步一次”的层；
- 状态注册表的间接访问 → view/span，kernel 内只读裸指针。

### 8.2 AMReX/GPU

- MultiFab 生命周期与 view 的边界 → 状态归 Core，块只拿 view；
- 新状态的 ghost/边界/AMR 同步 → 统一注册钩子。

### 8.3 迁移

- 现有 Evolve 顺序敏感 → 适配器阶段保留原顺序；
- `FieldType` 平铺枚举迁移 → 双写 + 编译期映射；
- 回归覆盖面 → 每阶段 P0/E3/S5/E6/L4/L5/L6。

### 8.4 物理正确性

- 任何新过程都要有守恒/熵/正性检查；
- 材料过程的热力学闭合与界面条件；
- 相互作用的自作用/边界/奇异性；
- 约束的投影是否破坏其他守恒。

### 8.5 过度设计

- 只为**已知角色**建接口；不为假想需求建抽象；
- 每个抽象配真实用例；允许“先用适配器、后抽接口”。

### 8.6 组合爆炸

- 插件之间的顺序/耦合由 PluginSpec 决定，但要限制“任意组合”的自由度；
- 提供**参考组合**与**守恒审计**，避免组合错误。

### 8.7 组织与流程

- 接口变更评审；
- 新物理 = 接口 + 注册 + 契约测试 + 文档；
- 迁移期保留“旧路径可回退”。

---

## 9. 参考

1. MOOSE 概览：<https://icmed.engin.umich.edu/wp-content/uploads/sites/176/2018/05/Summer-School-MOOSE-talk-2018-reduced.pdf>
2. Smilei：<https://arxiv.org/html/1702.05128v1>
3. PIConGPU：<https://github.com/computationalradiationphysics/picongpu>
4. Gkeyll 两流体（5/10 矩）：<https://ar5iv.labs.arxiv.org/html/1903.09618>
5. 电场诱导相分离（相场/介电混合）：<https://pubs.aip.org/aip/jcp/article/161/4/044704/3304666/Electric-field-induced-phase-separation-and>
6. C6F12O 沸腾/相变：<https://ouci.dntb.gov.ua/en/works/42B0kXV4/>
7. 数值宇宙学（PM 引力）：<https://indico.iap.fr/event/25/contributions/223/attachments/153/278/LesHouches_Lecture1.pdf>
8. WarpX 概览：<https://indico.slac.stanford.edu/event/8597/contributions/8188/attachments/3935/10601/2024_beam_beam.pdf>

---

## 10. 建议的决策清单（留给负责人，不在本文拍板）

| # | 决策点 | 选项 | 备注 |
|---|---|---|---|
| 1 | 状态注册表 | 保留 `FieldType` / typed registry / 双写过渡 | 影响语义与迁移成本 |
| 2 | 多态方式 | 静态多态优先 / 运行时虚函数 / 混合 | 性能 vs 扩展便利 |
| 3 | 统一等离子体块 | 粒子+流体同块 / 分块但共享接口 | 影响物种/EOS/摩擦复用 |
| 4 | 调度模型 | 固定顺序 / PluginSpec DAG / 手写 splitting | 组合能力 vs 可预测性 |
| 5 | 组合自由度 | 任意组合 / 受限组合 + 参考模板 | 组合爆炸风险 |
| 6 | 迁移范围 | 全量 / 只新物理走新接口 | 回归风险 |
| 7 | 1D 测试床与 3D | 先 1D 后 3D / 直接 3D | 验证速度 |
| 8 | 旧路径淘汰时机 | 立即 / 迁移完成后 | 维护成本 |
| 9 | 性能预算 | 抽象层开销上限 / 热点白名单 | 需基准 |
| 10 | 契约测试范围 | 每插件 / 每角色 / 抽样 | 测试成本 |

---

## 11. 一句话

**建议把“新物理”归入有限个通用角色（场方程、相互作用、物种模型、材料模型、约束、源、边界、诊断、积分器），用统一接口 + PluginSpec 声明 + 工厂注册 + typed 状态 + 契约测试来支撑任意组合；引力和相变只是“单角色 1–2 个文件”的成本示意，不是架构针对的对象；迁移用适配器 + strangler + 逐位回归；所有取舍列在 §10，交由负责人决定，本文不拍板。**


---

## 12. 自检与修正（对本文建议的审查）

> 以下是自查发现的问题与修正；前文相关段落若与本節冲突，以本節为准。

### 12.1 已修正的结构性错误

1. **“DAG 调度”不成立（最严重）**
   - 场与源是**循环依赖**（场依赖源、源依赖场）；reads/writes 的 DAG 无法给循环排序。
   - 修正：DAG 只用于**同一 split 阶段内的偏序**；循环由 operator split / IMEX / 定点迭代显式处理（§5.3 已改）。

2. **`IInteraction::solve` 与 `IFieldEquation` 重复**
   - 相互作用不应自己解场；它只负责 deposit/gather，并声明使用哪个场方程。
   - 修正：`IInteraction` 去掉 `solve`，加 `field_equation()`（§5.2 已改）。

3. **`IMaterialModel` 三合一破坏模块性**
   - 物性（ε/σ，场求解器要）、状态过程（相变/反应）、响应（σE/ρ_medium，源要）的消费者不同；
     合并会强迫材料模型实现全部方法。
   - 修正：拆回 `IMaterialProperty` / `IMaterialProcess` / `IMaterialResponse`（§5.2 已改）。

4. **`ISpeciesModel` 的 `push`/`transport` 语义重叠**
   - PIC 的粒子推与流体的输运是同一语义；`push` 与 `transport` 并存会让人误解。
   - 修正：统一为 `advance`（§5.2 已改）。

5. **“open slot”说法过强**
   - 编译期 `FieldKey<R,T>` 只能覆盖已知状态；插件自定义的**运行时命名 slot** 需要类型擦除，
     二者要分清（§5.4 已注）。

6. **“1–2 个 .cpp”低估**
   - 除核心实现外，还有状态 slot、边界、输入、诊断、契约测试等胶水；
   - 真实成本应表述为“**1–2 个核心文件 + 胶水**”（§6 已改）。

7. **“任意组合”不能靠声明保证守恒/一致**
   - `conservative` 标志只是声明，不是证明；组合可能重复计数或破坏守恒；
   - 修正：必须配守恒审计、参考组合、契约测试；不承诺“自动正确”（§8.6 已提，此处强调）。

8. **性能结论偏乐观**
   - “静态多态 + 虚函数每步一次”只解决**每步**的分派；每个 cell 的插件代码仍需
     **templated kernel + 类型擦除的 kernel launcher**，并考虑 kernel launch/融合/内存布局；
     需要基准，不能只靠声明。

9. **“P0 逐位不变”是目标而非保证**
   - 适配器保持调用顺序通常逐位；但新接口的浮点顺序可能变化；
   - 修正：以“逐位优先、数值回归为兜底”，并记录偏差。

10. **文献引用**
    - MOOSE/Smilei/PIConGPU/Gkeyll 的“模式”描述是概括，不是逐条原文；
    - 相变/PM 引力的具体公式与假设需按原文核对后再用于实现。

### 12.2 仍然成立的部分

- 现状审计（上帝类、平铺 FieldType、介质嵌在 FDTD、两套流体、Evolve 写死顺序）；
- 通用角色分类（场/耦合/物种/材料/约束/源/边界/诊断/积分器）；
- 插件 + 注册 + 输入 + 契约测试的总体方向；
- 适配器 + strangler + 分阶段迁移；
- 风险清单（性能、AMReX、迁移、物理正确性、过度设计、组合爆炸）。

### 12.3 需要补充的决策

- 循环耦合的**具体 splitting/IMEX 方案**（Strang 顺序、场源冻结、迭代次数、收敛判据）；
- **守恒审计**的机制（能量/动量/电荷/熵）与容差；
- 运行时命名 slot 的**类型擦除**实现（string/id → 布局/访问器）；
- **kernel launcher** 的设计（templated kernel + 函数指针/变体）；
- 组合的**受限集合**与参考模板（避免任意组合的不可控）。

---

## 13. 一句话（修正版）

**通用角色 + 插件注册 + typed/命名状态 + 契约测试的方向成立；但“DAG 调度”“interaction 自带 solve”“材料三合一”“1–2 个 .cpp 且不碰核心”“任意组合自动正确”这几条需要修正为：循环依赖用 split/IMEX 显式处理、interaction 只做源/力并引用场方程、材料拆成物性/过程/响应、核心文件之外有胶水、组合必须配守恒审计与受限模板。**

---

## 14. 审阅结论的采纳（对 `ARCHITECTURE_REVIEW.md`）

> 审阅文档：`Docs/research_notes/ARCHITECTURE_REVIEW.md`（提交 `aef3175`）。
> 以下三条**全部采纳**；§5.4/§6/§5.2 的相关表述以本节为准。

### 14.1 三处修正

1. **typed registry 不是必经之路（采纳）**
   - 事实：`ablastr::fields::MultiFabRegister` 已经是**字符串键**的注册表；新状态今天就能注册，`FieldType` 只是已知字段的便捷枚举。
   - 修正：**保留 `FieldType` + 复用字符串注册表**；不新造 typed registry；等出现“同名冲突/拼写错误”的真实痛点再考虑。
   - 影响：§5.4 / §12.5 的“typed key”从“必经”降为“可选、后置”。

2. **“P0 逐位不变”做不到（采纳）**
   - 事实：`EvolveB` 的 `SubcyclingHalf::FirstHalf` 半格、PSATD 的 Vay J 与 rho 半格对齐、适配器改变 `fillBoundary` 时机会改变浮点。
   - 修正：P0 的验收改为**数值门（相对误差 ≤1e-12）+ 少数简单路径逐位**；基线用 `P0 0.0927/0.8906`、`S5 0.619830`、`E6 0.100042`。
   - 影响：§6 P0、§12.9 的“逐位优先”改为“数值门为准，逐位为锦上添花”。

3. **接口粒度与 GPU 形状（采纳并补形状）**
   - 事实：`solve(FieldState&, dt)` 是宿主侧抽象；插件上 GPU 需要 `Array4`；且未体现 E/B 分开推进 + 半格交错。
   - 修正（补进设计）：

```cpp
struct FieldView { amrex::Array4<amrex::Real> a; amrex::Box bx; int ncomp; };

class IFieldEquation {
public:
    virtual void advance(FieldState&, const FieldSource&, Real dt, TimeLayer layer) = 0;
    // layer ∈ {Full, FirstHalf, SecondHalf}
    virtual void boundary(FieldState&) = 0;
    virtual ~IFieldEquation() = default;
};
// 插件实现：template<class Impl> launch(FieldView, ...)（内联、无虚函数）
// 虚函数只在“选择 launcher”这一层（每步一次）
```

   - 影响：§5.2 的 `solve` 改为 `advance(..., TimeLayer)`；§12.8 的 launcher 补上 `FieldView` 形状。

### 14.2 采纳的 MVP 路线（medium 单块，2–3 周）

```
只做 medium 块接口化：
  IMediumProperty  —— 给场求解器 fill ε/σ/μ（适配现有 MacroscopicProperties）
  IMediumSource    —— 给源算 ρ_medium / σE（适配现有 ⟨J_c⟩ 路径）
  回归门：L0–L6（P0/S5/E6 基线；数值门 ≤1e-12）
不动 FieldType / Evolve / 场求解器（适配器 + 双写）
等第二个真实需求（μ 张量 / 第三材料）再决定是否推广
```

理由（采纳审阅意见）：

- 唯一“有真实用例、有现成回归、不阻塞主线”的切片；
- 验证 §12.1(3) 材料三分法是否够用——medium 的真实结构正是 `ε/σ` 给场、`ρ_medium` 给源、界面给 BC；
- 成本 2–3 周，远低于全量 12–21 人月。

### 14.3 对 §10 决策清单的采纳值

| # | 决策点 | 采纳 |
|---|---|---|
| 1 | 状态注册表 | 保留 `FieldType` + 复用字符串注册表（不新造 typed registry） |
| 2 | 多态方式 | 静态多态优先；虚函数只在 launcher 选择层 |
| 3 | 统一等离子体块 | 分块但共享接口（不先合并） |
| 4 | 调度模型 | 固定顺序 + 显式 split 阶段（不用 DAG） |
| 5 | 组合自由度 | 受限组合 + 参考模板 |
| 6 | 迁移范围 | 只新物理走新接口（medium MVP） |
| 7 | 1D/3D | 先 1D，后 3D |
| 8 | 旧路径淘汰 | 第二个真实需求出现后再议 |
| 9 | 性能预算 | 先基准，后定上限 |
| 10 | 契约测试 | 每插件 + 每角色抽样 |

### 14.4 一句话（采纳审阅后）

**方向成立，但全量 P0–P5 不划算；采纳审阅：保留字符串注册表、P0 用数值门而非逐位、接口补 `TimeLayer` + `FieldView` 形状；先做 medium 单块 MVP（`IMediumProperty` + `IMediumSource`，2–3 周，L0–L6 回归），不动核心，等第二个真实需求再推广。**
