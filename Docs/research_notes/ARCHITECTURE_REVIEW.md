# 对《ARCHITECTURE_MODULARIZATION_PROPOSAL.md》的可行性与正确性审阅

> 审阅对象：`Docs/research_notes/ARCHITECTURE_MODULARIZATION_PROPOSAL.md`（498 行，含 §12 自检）
> 方法：逐条对照 `Source/` 实际代码；凡本文给出的数字均来自 `grep/wc` 实测。
> 结论分三档：**成立 / 需修正 / 不成立**。

---

## 0. 结论摘要

| 维度 | 评价 |
|---|---|
| **方向** | 成立。"通用角色 + 插件注册 + 接口/适配器 + strangler"是主流做法，与 MOOSE/Smilei/Gkeyll 的实践一致 |
| **§3 现状审计** | 基本成立（4/5 条与代码一致，1 条数字略偏） |
| **§12 自我修正** | 4 条主要修正**全部正确**，是本稿最有价值的部分 |
| **可行性（全量 P0–P5）** | **低**：600+ 源文件、5 条场求解路径，P0 一项保守 3–6 人月，且"逐位不变"目标不现实 |
| **可行性（medium 单块 MVP）** | **高**：2–3 周可落地，且有现成回归（L0–L6）作安全网 |
| **是否建议现在做** | **不建议全量**。建议只对 medium 块做接口化试点，等第二个真实需求出现再推广 |

---

## 1. 事实核对（对照代码）

| 文档论断 | 实测 | 判定 |
|---|---|---|
| `FieldType` 平铺枚举 ~100 项 | **76** 项（`Source/Fields.H`，`AMREX_ENUM(FieldType, ...)`） | ✅ 成立（数字略高） |
| `WarpX` 是上帝类 | `WarpX.H` **1808 行**、约 **305** 个成员声明 | ✅ 成立 |
| `MacroscopicProperties` 嵌在 FDTD 目录 | 实际在 `Source/FieldSolver/FiniteDifferenceSolver/MacroscopicProperties/` | ✅ 成立 |
| 两套流体不共享 EOS/闭合/摩擦 | `Source/Fluids/` 内 **0** 处引用 `EquationOfState`/`FluidFriction`；`MomentsSolver/Coupling/` 独有 `TenMomentClosure.H`、`MomentCoupling.H` | ✅ 成立 |
| `Evolve()` 物理顺序写死 | `WarpXEvolve.cpp` 内 `EvolveF(0.5dt)` → `EvolveB(0.5dt, SubcyclingHalf::FirstHalf)` → `EvolveE` / `MacroscopicEvolveE` 分支写死 | ✅ 成立 |
| `MomentsSolver/Equations` 是 equation-object 雏形 | `Euler5Eqn.H`、`TenMomentEqn.H`、`MomentEquation.H`、`MaxwellEqn.H` | ✅ 成立 |

**结论：§3 的现状审计可信**，可作为设计依据。

---

## 2. 正确性：文档遗漏/高估的三处

### 2.1 【高估】"typed registry 是必经之路"——**不成立**

文档 §5.4 / §12.5 认为新状态需要"typed key + 类型擦除的 registry"，并把它列为待决策项。

**代码事实**：`ablastr::fields::MultiFabRegister` **已经是字符串键**：

```cpp
// MultiFabRegister.H
/** This registers a new MultiFab under a unique name, allocates it and ... */
```

即：**今天就能在不改 `FieldType` 的前提下注册任意新状态**（`register` + 名字 + 布局）。`FieldType` 只是**已知字段的便捷枚举**，不是注册表的硬限制。

- 影响：文档把"状态注册"低估为需要新造的组件，实际上**先做适配器层即可复用**；
- 建议：`FieldType` 保持不动，新插件直接用字符串名注册；等出现"同名冲突/拼写错误"的真实痛点再考虑 typed key。

### 2.2 【不成立】"P0 阶段现有 regression 逐位不变"——**做不到**

文档 §7 P0 写"验收：现有 regression 逐位/数值不变"，§12.9 已部分修正为"逐位优先、数值兜底"，但仍偏乐观。

原因（代码层面）：

```
1. EvolveB 带 SubcyclingHalf::FirstHalf（半格 B），接口若写成 solve(FieldState, dt) 会掩盖时间层；
2. PSATD 路径里 J 由 Vay 构造、且与 rho 的半格对齐耦合（WarpXEvolve.cpp 注释明说）；
3. 适配器只要改变一次 fillBoundary 的时机或顺序，浮点结果就变。
```

- 正确表述：**数值一致（相对误差 ≤1e-12）+ 少数简单路径逐位**；逐位只能作为"锦上添花"，不能作为 P0 的验收门。
- 我们已有的基线正好可用作数值门：`P0 0.0927/0.8906`、`S5 0.619830`、`E6 0.100042`。

### 2.3 【缺口】接口粒度与 GPU kernel 的冲突——**未给出可落地的形状**

§5.2 的接口是：

```cpp
virtual void solve(FieldState&, const FieldSource&, Real dt) = 0;
```

问题：`FieldState&` 是**宿主侧抽象**。插件要在 GPU 上跑，必须自己拿到 `Array4`；而 `FieldState` 若隐藏了 `MultiFab`，插件就得回调 Core 取数组——这会让"每步一次虚函数"退化成"每 kernel 一次虚函数 + 类型擦除"，正是 §12.8 担心的，但文档没给出解法。

**可落地的形状**（建议补进文档）：

```cpp
// 1) Core 提供“按需取裸视图”的窄接口（不暴露 MultiFab 语义）
struct FieldView { amrex::Array4<amrex::Real> a; amrex::Box bx; int ncomp; };
// 2) 插件实现是模板/静态多态，注册时绑定到一个 kernel launcher
template <class Impl> void launch(FieldView, ...);   // 内联，无虚函数
// 3) 虚函数只出现在 “选择 launcher” 这一层（每步一次）
```

另外一条硬约束：`IFieldEquation::solve` 与现有 **E/B 分开推进 + 半格交错** 的结构不匹配。建议改成：

```cpp
virtual void advance(FieldState&, const FieldSource&, Real dt, TimeLayer layer) = 0;
// layer ∈ {Full, FirstHalf, SecondHalf}
```

否则 Maxwell 的时间层语义会被藏进实现，HybridPIC/PSATD 无法接。

---

## 3. 文档中**正确且有价值**的部分

1. **§12.1(1) DAG 调度不成立**：场与源循环依赖，reads/writes 的 DAG 无法排序。这是原稿最严重的错误，修正正确。
2. **§12.1(2) `IInteraction` 去掉 `solve`**：相互作用只做 deposit/gather + 声明所用场方程，正确。
3. **§12.1(3) 材料拆成 Property/Process/Response**：三者消费者不同（场求解器 / 状态演化 / 源），合并会强迫实现全部方法。**这一条与我们 medium 模块的实际结构完全吻合**（`ε/σ` 给场、`ρ_medium` 给源、界面给 BC），是本文最有实操价值的一条。
4. **§12.1(7) `conservative` 只是声明不是证明**：正确，必须配守恒审计。
5. **§8 风险清单**（性能、AMReX 生命周期、迁移、物理正确性、过度设计、组合爆炸）覆盖完整。

---

## 4. 可行性评估

### 4.1 规模与成本

```
Source/：626 个源文件（368 .H + 258 .cpp）
场求解路径：FDTD(Yee/CKC) · PSATD · 静电 · 隐式 · HybridPIC  —— 5 条
```

| 阶段 | 内容 | 保守估计 | 主要风险 |
|---|---|---|---|
| P0 | 接口 + 适配器 + 调度器 | **3–6 人月** | 5 条路径都要接；时间层语义（2.3） |
| P1 | Plasma 块（EOS/闭合/摩擦） | 2–4 人月 | 两套流体统一会改 `Fluids/` 与 `MomentsSolver/` |
| P2 | Fields 块 + typed registry | 2–3 人月 | 与 2.1 重复投入（注册表已可用） |
| P3 | Medium 块 | **0.5–1 人月** | 已有 L0–L6 回归 |
| P4 | Coupling 层（split/IMEX + 守恒审计） | 3–5 人月 | 守恒审计是新功能 |
| P5 | 淘汰旧路径 | 1–2 人月 | 回归面最大 |

**全量约 12–21 人月**，且期间主线开发几乎停滞。

### 4.2 与当前项目需求的匹配度

当前实际需求集中在 **medium 块**（L0–L7 验证、σ 张量、界面电荷）。而 P0/P1/P2/P4 对当前需求**零收益**，只有 P3 直接相关。

---

## 5. 建议（MVP 路线）

```
1. 只做 medium 块接口化（P3 的缩小版），2–3 周：
   - IMediumProperty：给场求解器 fill ε/σ/μ（对接现有 MacroscopicProperties）
   - IMediumSource  ：给源算 ρ_medium / σE（对接现有 ⟨J_c⟩ 路径）
   - 现有 L0–L6 回归作数值门（P0/S5/E6 基线）
2. 不动 FieldType、不动 Evolve、不动场求解器：适配器 + 双写
3. 等第二个真实需求（例如 μ 张量、或第三材料）出现，再决定是否推广到 Fields/Plasma
```

**理由**：这是唯一"有真实用例、有现成回归、且不阻塞主线"的切片；同时它验证了文档 §12.1(3) 的材料三分法是否真的够用。

---

## 6. 对 §10 决策清单的建议值

| # | 决策点 | 建议 |
|---|---|---|
| 1 | 状态注册表 | **保留 `FieldType` + 复用字符串注册表**（不新造 typed registry） |
| 2 | 多态方式 | 静态多态优先，虚函数只在 launcher 选择层 |
| 3 | 统一等离子体块 | **分块但共享接口**（不要先合并） |
| 4 | 调度模型 | 固定顺序 + 显式 split 阶段；**不用 DAG** |
| 5 | 组合自由度 | **受限组合 + 参考模板** |
| 6 | 迁移范围 | **只新物理走新接口**（strangler） |
| 7 | 1D 测试床与 3D | 先 1D（`MomentsSolver` 已是） |
| 8 | 旧路径淘汰 | 迁移完成后 |
| 9 | 性能预算 | 需基准；暂定抽象层 ≤2% |
| 10 | 契约测试范围 | 每插件 |

---

## 7. 一句话

**方向成立、§3 审计可信、§12 的四条修正（DAG/interaction/材料三分/语义统一）都对；但"typed registry 必经"被高估（注册表已是字符串键）、"P0 逐位不变"不现实（半格 B / PSATD 的 J 对齐）、接口粒度与 GPU 的形状缺失；全量 P0–P5 约 12–21 人月且与当前需求不匹配，建议只做 medium 块的接口化 MVP（2–3 周，用 L0–L6 作数值门），其余等第二个真实需求再议。**
