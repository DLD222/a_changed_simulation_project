# 只增不减的模块插入架构：调研与建议（建议稿，非决定）

> 需求：**旧源码只增不减、只插入模块**——不改旧文件已有逻辑，不删旧路径；
> 新物理以“模块”插入。本文调研可行的插入方式，给出选项、代价与风险；**不拍板**。
> 日期：2026-11-04 · 基于 `~/WarpX_fluid` 源码与公开资料。

---

## 0. 先明确“只增不减”的两种含义

| 级别 | 含义 | 旧文件改动 |
|---|---|---|
| **A. 零改旧文件** | 旧 `.cpp/.H` 一行不改；只新增文件/新目标/新目录 | 0 |
| **B. 只做加法** | 旧文件只**新增**行（新 case、新调用、新注册），已有行不动、不删 | 少量新增 |

本文把 A 作为首选目标，B 作为“旧特性不能丢”时的备选；两者都**不删、不改旧逻辑**。

---

## 1. 代码事实：旧 WarpX 已经暴露了“可插入”的公共步骤

`Source/WarpX.H`：

```
public:                       // line 89
  void InitData();            // line 131
  void Evolve(int numsteps);  // line 133
  void EvolveE(dt, start_time);                 // 505
  void EvolveB(dt, SubcyclingHalf, start_time); // 507
  void EvolveF/EvolveG(...);                    // 509–516
  void MacroscopicEvolveE(...);                 // 518–520
  void PushParticlesandDeposit(...);            // 694/700
protected:                    // line 1128
```

这意味着：**新的驱动（新 main/新 app）可以调用 `InitData` + 这些公共步骤，自己写时间循环，并在步骤之间插入新模块，而完全不改旧文件**。

同时，旧代码的“工厂”都是显式 switch/if：

```
MultiParticleContainer: PhysicalParticleContainer / RigidInjected / Photon / Laser
MultiFluidContainer   : WarpXFluidContainer
FieldSolver           : FiniteDifferenceSolver / SpectralSolver
MultiDiagnostics      : Full / BT / BoundaryScraping
```

新增类型需要在这些工厂里加一个 **case**（B 级“只增”），否则就要走 A 级的新驱动。

---

## 2. 文献与工程模式

| 模式 | 机制 | 参考 |
|---|---|---|
| **微内核/插件** | 核心最小 + 插件注册 + 扩展点；插件通过 SPI 接入 | [Microkernel (Plugin) Architecture](https://github.com/archman-dev/website/blob/main/docs/architectural-styles/microkernel-plugin.mdx) |
| **MOOSE App** | 框架不动；新物理写成新 App + 新 Kernel，注册进 App | [MOOSE 自定义 Kernel](https://raw.githubusercontent.com/idaholab/moose/2b29e9a5ed7a4e8e5811c3acc927e9e82cdcd91f/modules/doc/content/examples/ex02_kernel.md)、[MOOSE](https://en.wikipedia.org/wiki/MOOSE_(software)) |
| **动态插件** | 插件编译为 `.so`，运行时 `dlopen` 加载；注册表在插件/宿主 | [dlopen 依赖识别](https://lwn.net/Articles/970803/) |
| **静态注册** | 插件源文件的静态构造器把工厂注册进全局表；链接时 `--whole-archive` 保证不被丢弃 | 通用 C++ 插件模式 |
| **WarpX 公共步骤** | `InitData`/`EvolveB`/`EvolveE`/`PushParticlesandDeposit` 均为 public | 本项目源码 |

---

## 3. 可行方案（按旧码改动量排序）

### 方案 A：**Sidecar 新驱动（零改旧文件）**——首选

```
新增：
  Source/Multiphysics/            # 新层：接口、注册、调度、状态视图、插件
  Apps/Multiphysics/main.cpp      # 新可执行：InitData() + 自写循环 + 插入新模块
  CMake/Make.package 新增目标     # 只加新 target/目录
旧码：
  0 行改动；旧 main 保留；旧 WarpX 类原样
```

新驱动的循环 = 旧 `Evolve`/`OneStep_nosub` 的**公共步骤**复制 + 新模块插入：

```
InitData()
loop:
  PushParticlesandDeposit(...)
  EvolveB(0.5dt, FirstHalf)
  EvolveF/G, PML, FillBoundary
  [新模块: Multiphysics::step()]
  EvolveB(0.5dt, SecondHalf)
  EvolveE / MacroscopicEvolveE
  diagnostics / sync
```

- 优点：**旧源码零改动**；新物理完全在新层；旧 `main`/旧算例不受影响；
- 代价：新驱动要复制旧循环的步骤（含 PML/EB/AMR/诊断/subcycling 的处理）；旧循环以后的变化不会自动同步；
- 适用：以新物理为主、旧特性可接受“复制版”的场合。

### 方案 B：**只做加法的扩展点（旧文件只加新行）**

```
旧文件新增（不删不改）：
  MultiParticleContainer/ MultiFluidContainer/ FieldSolver/ MultiDiagnostics:
    加一个 case/注册（新类型）
  Evolve/WarpXEvolve.cpp:
    在现有步骤之间新增一次 Multiphysics::step()（默认 Off）
  新文件：Source/Multiphysics/**（接口、注册、插件）
```

- 优点：旧循环/旧特性/旧算例全部保留；新模块注册即可用；
- 代价：旧文件有**新增行**（B 级“只增”）；默认 Off 保证旧路径不变；
- 适用：旧特性（AMR/EB/PML/Python/诊断）不可丢、且能接受少量新增行的场合。

### 方案 C：**动态插件（`.so` + dlopen）**

```
新模块编译为 libmultiphysics_plugin.so；新驱动（方案 A）或宿主加载；
注册表在新层；旧码 0 改动（若宿主是新驱动）
```

- 优点：可运行时插拔；旧码 0 改动（配合方案 A）；
- 代价：ABI/版本/符号可见性；GPU/多进程下的加载；调试复杂；
- 适用：需要“不改二进制就换模块”的场合；一般项目不必先上。

### 方案 D：**Python 层扩展**

```
pywarpx 层新增 Python 模块，调用旧 C++ 的公开接口；旧 C++ 0 改动
```

- 优点：最快验证；旧码 0 改动；
- 代价：性能/GPU 受限；只适合原型/分析，不适合生产内核。

---

## 4. 对比

| 维度 | A（Sidecar） | B（只加新行） | C（动态插件） | D（Python） |
|---|---|---|---|---|
| 旧源码改动 | **0** | 少量新增行 | 0（配合 A） | 0 |
| 旧特性保留 | 复制版（需同步） | ✅ 原样 | 复制版 | 复制版 |
| 新物理扩展 | 新层插件 | 新层插件 + case | 插件 | Python |
| 性能 | 原生 | 原生 | 原生 | 受限 |
| 复杂度 | 中（复制循环） | 低 | 高 | 低 |
| 适用 | 新物理为主 | 旧特性不可丢 | 运行时可插拔 | 原型 |

---

## 5. 建议路线（建议，非决定）

```
S0（2–3 周）：方案 A 骨架
  - Apps/Multiphysics/main.cpp：InitData + 自写循环（先 1D/2D，简单边界）
  - Source/Multiphysics：接口 + 注册 + 调度 + FieldView/TimeLayer
  - 第一个模块：medium（ε/σ/μ + ρ_medium，适配 MacroscopicProperties）
  - 回归：L0–L6 数值门（P0/S5/E6 基线）
旧码：0 改动
```

```
S1：若旧特性（AMR/EB/PML/诊断）必须保留，切换/叠加方案 B：
  只在工厂/Evolve 加新行（默认 Off）；旧路径不变
```

```
S2：第二、第三个新模块（EOS/摩擦、任意 interaction/field）——验证“只增”扩展
S3：需要运行时可插拔时，再评估方案 C
```

---

## 6. 风险与边界

```
- 方案 A 的循环复制：旧循环的 PML/EB/AMR/subcycling/诊断必须逐项对齐，否则数值/功能缺失
- 旧码演进：旧 Evolve 改动不会自动进新驱动；需要同步机制（测试/对比）
- 方案 B 的“新增行”：必须默认 Off，且旧路径逐位/数值不变
- 动态插件：ABI/符号/GPU/多进程；先不做
- 测试：旧路径回归 + 新模块契约测试 + 组合守恒审计
- 本文件只是建议；A/B 的取舍由负责人决定
```

---

## 7. 参考

1. Microkernel (Plugin) Architecture：<https://github.com/archman-dev/website/blob/main/docs/architectural-styles/microkernel-plugin.mdx>
2. MOOSE 自定义 Kernel：<https://raw.githubusercontent.com/idaholab/moose/2b29e9a5ed7a4e8e5811c3acc927e9e82cdcd91f/modules/doc/content/examples/ex02_kernel.md>
3. MOOSE：<https://en.wikipedia.org/wiki/MOOSE_(software)>
4. dlopen 依赖识别：<https://lwn.net/Articles/970803/>
5. WarpX 概览：<https://indico.slac.stanford.edu/event/8597/contributions/8188/attachments/3935/10601/2024_beam_beam.pdf>
6. 本项目源码：`Source/WarpX.H`（public 步骤）、`Source/Particles/MultiParticleContainer.cpp`、`Source/Fluids/MultiFluidContainer.cpp`、`Source/Diagnostics/MultiDiagnostics.cpp`

---

## 8. 待决策清单（留给负责人）

| # | 决策点 | 选项 | 备注 |
|---|---|---|---|
| 1 | 旧码改动级别 | A 零改动 / B 只加新行 | 本文建议 A 起步 |
| 2 | 新驱动范围 | 1D/2D 最小循环 / 完整复制旧循环 | 复制越全，旧特性越全 |
| 3 | 旧特性优先级 | AMR/EB/PML/诊断哪些必须先有 | 决定 A 还是 B |
| 4 | 动态插件 | 现在 / 以后 | 一般不必先上 |
| 5 | 新层位置 | `Source/Multiphysics` / `Apps` | 影响构建与依赖 |
| 6 | 测试策略 | 旧回归 + 新契约 / 只新测试 | 建议前者 |

---

## 9. 一句话

**旧 `WarpX` 的 `InitData`/`EvolveB`/`EvolveE`/`PushParticlesandDeposit` 都是 public，因此“零改旧文件”的 sidecar 新驱动可行：新增 `Apps/Multiphysics` + `Source/Multiphysics`，用旧公共步骤自写循环并插入新模块；若旧特性（AMR/EB/PML/诊断）不能丢，则退到“只加新行”的方案 B（工厂/Evolve 各加一个默认 Off 的入口）；动态插件与 Python 层作为后续/原型选项。本文只给选项与代价，不拍板。**
