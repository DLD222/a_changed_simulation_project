# 电磁 / 介质 / EB 验证阶梯：10 个算例

> 面向**磁流体（矩方程）仿真**的能力验收、误差标定与回归计划
>
> - 适用代码：`<repo>`（WarpX 26.10 分支 + `Source/MomentsSolver`）
> - 版本：v2.0（草案；C2 已按「与 C1 同物理、仅几何实现不同」重写）　|　日期：2026-10-31
> - v2 变更：C2 由「PEC cut-cell」改为「与 C1 完全相同的介质算例，仅界面由 cut-cell + 亚网格平均表示」；
>   PEC 版本降级为 **C2-0 极限自检**（现在即可跑）。M4（AMR）的实施调研见
>   `M4_AMR_IMPLEMENTATION_RESEARCH_v1.md`。
> - **状态声明**：本文只做设计与规划，**尚未运行任何算例**。文中出现的数字只有两类：
>   ① 解析解/理论值；② **预注册**的验收阈值。**没有任何仿真结果数字**。
> - 本文档是"先查论文再跑、严禁数据造假"这两条规则的操作化版本：每一节都给出文献门槛、
>   解析解来源、数据来源与留痕方式。

---

## 0. 一页速览

### 0.1 为什么磁流体代码要先做电磁/介质/EB 验证

本项目要开发的是"动理学离子 + 矩（流体）电子 + 全麦克斯韦"路线。矩方程本身（5 矩 / 10 矩）
只描述物质，真正把物质和场耦合起来的是麦克斯韦求解器、材料（ε、μ、σ）本构、以及
几何（以 cut-cell EB 表示的介质/电极）。**场求解器错 1%，等离子体频率和波模就会错 1% 以上**，
而这类错误不会在"看起来像"的场图里暴露。因此必须先把场/材料/几何这条链在**有解析解的
标准问题上**标定到已知误差，再谈 MHD。

本阶梯从"1D 平面波打平面界面"走到"核壳球 + 曲面 cut-cell + AMR"，每一级都只引入一个新自由度
（角度、曲率、材料数、损耗、多层、倏逝、等离子体、磁化、自适应）。

### 0.2 十个算例一览

| # | 算例 | 核心理论解 | 主要验证对象 | 前置能力 |
|---|------|-----------|-------------|---------|
| C1 | 平面介质界面：Fresnel + Brewster | 复 Fresnel 公式、Snell、$\tan\theta_B=n_2/n_1$ | 材料求解器、斜入射几何 | 现成（3D/2D/1D 构建） |
| C2 | 斜平面介质界面 cut-cell（**与 C1 同物理，仅几何实现不同**） | 复 Fresnel（同 C1）；PEC 极限 $r=-1$（C2-0） | cut-cell 几何、亚网格平均、几何误差 | C2-0 现成（ECT）；C2 需 M2 |
| C3 | 介质平板 Fabry–Perot | Airy 公式、FSR、finesse | 薄层、多次反射、能量守恒 | 现成（材料） |
| C4 | 圆柱与球的 Mie 散射 | Mie 级数、柱 Bessel/Hankel 级数 | 曲面几何、散射截面、近场 | PEC 现成；介质需 M2 |
| C5 | 有耗介质与趋肤深度 | 复折射率、$\delta_s=\sqrt{2/(\omega\mu\sigma)}$ | 损耗项 α/β、吸收率 | 现成（σ） |
| C6 | 多层膜：增透膜与 Bragg 镜 | 特征矩阵（传输矩阵） | 亚网格薄层、多界面 | 厚层现成；薄层需 M2 |
| C7 | 全内反射与倏逝波 | $\sin\theta_c=n_2/n_1$、衰减长度 $1/\kappa$ | 倏逝场、PML 吸收 | 现成（材料 + PML） |
| C8 | 无磁化等离子体色散 | $\omega^2=\omega_p^2+c^2k^2$、截止、趋肤 | 电荷/电流沉积、色散测量 | 现成（PIC / 冷流体） |
| C9 | 磁化等离子体波模 | Stix 色散、whistler、Alfvén | 混合/Darwin 求解器、磁化物理 | 现成（hybrid / Darwin） |
| C10 | 核壳球（Aden–Kerker）+ 曲面 cut-cell + AMR | 双层 Mie 级数 + 两个极限 | 综合：曲率 + 多材料 + AMR + 守恒 | 需 M2（介质 EB）+ M4（AMR 材料） |

“M2/M4”是既有的分阶段计划编号（见 §8.3）：M2=电介质 EB 单层，M4=AMR。

### 0.3 三条硬规则（违反任一条，算例作废）

1. **先文献后计算（G0）**：每个算例在跑第一次之前，必须提交"文献卡"——列出所用公式的
   出处（含公式编号/页码/DOI）、约定（时间因子 $e^{-i\omega t}$ 还是 $e^{+i\omega t}$、
   $\varepsilon$ 的虚部符号、$\mu_r$ 是否为 1）、以及**要复现的文献数值表**。没有文献卡不许跑。
2. **一个数一个来源（G2）**：报告里出现的每个数字，必须能由 `analyze.py` 从
   `raw/` + `theory.py` 重新算出来，并给出文件哈希。手工改数、口头传数一律视为造假。
3. **失败算例不删**：不收敛、超差、发散的算例原样保留，写进 `STATUS.md`，并在报告里
   单列"未通过项"。删数据比错数据更糟。

### 0.4 与旧工作的关系（重要）

2026-10-30 的 `bench3d`/`FINDINGS.md` 工作（26 例斜入射、EB 介质五轮）**其源码树与全部数值
已在 2026-10-30 23:57 之后被删除**（磁盘现状已核查：`WarpX_changed` 为 2026-10-31 00:06 的
全新 clone，旧树与 `bench3d/` 不存在）。因此：

- 旧报告里的数字**不可复查**，本计划**不引用、不复用**任何旧数值作为基线；
- 旧报告在方法论上仍然有效的三条经验被吸收进本计划：
  (i) 60° 异常来自单平面 TF/SF 面源的受迫 $k_z\approx0$ 分量 → 本计划用**倾斜界面 + 法向源**
  代替斜源，从源头避免该假场；
  (ii) 真空对照（仪器地板）必须每个算例都做；
  (iii) 未填充域会让"能量增长"看起来像不稳定 → 判据用**时间窗内最大值**，不用全局能量。

---

## 1. 现状：这台机器上的代码到底能做什么（实测）

> 本节每条结论都可复核；给出文件与行号（行号以当前工作树为准）。

### 1.1 求解器与维度

| 项 | 实测 | 证据 |
|---|---|---|
| 电磁求解器 | `yee`（默认）、`ckc`、`psatd`、`ect`、`hybrid`（Darwin/Ohm） | `Source/Utils/WarpXAlgorithmSelection.H`（`AMREX_ENUM(ElectromagneticSolverAlgo, None, Yee, CKC, PSATD, ECT, HybridPIC, ...)`） |
| 宏观材料算法 | `BackwardEuler`（默认）、`LaxWendroff` | 同上（`MacroscopicSolverAlgo`） |
| 维度 | 1D/2D/3D/RZ 均可构建；当前构建只配了 **3D** | `build/CMakeCache.txt`：`WarpX_DIMS:STRING=3` |
| 当前构建开关 | EB=ON，OMP，MPI=ON，PYTHON=**OFF**，openPMD=**OFF**，FFT=**OFF**，QED=OFF | `build/CMakeCache.txt` |

### 1.2 材料（宏观介质）

- `algo.em_solver_medium = macroscopic` 时启用；ε、μ、σ 可用常数或 `*_function(x,y,z)` parser
  （`Docs/source/usage/parameters.rst:3888`、`:4071-4101`）。
- **只允许单层**：`lev==0`，否则 abort（`Source/WarpX.cpp:2631-2632`：
  “Macroscopic properties are not supported with mesh refinement.”）。
- **RZ 禁用**：`Source/FieldSolver/FiniteDifferenceSolver/MacroscopicEvolveE.cpp:53`
  (“currently macro E-push does not work for RZ”)。
- 代码对 Cartesian 1D/2D/3D 编译（外层 `#if !defined(WARPX_DIM_RZ)`），但源文件注释写
  “currently works only for 3D”；**1D/2D 是否真的可用需在 G2 中实测确认**（本计划把这条列为
  C1 的第一个子任务）。
- ε、μ、σ 存于 cell center，再插值到 E/B 交错位置（`MacroscopicEvolveE.cpp:183-228`）；
  **没有任何 EB 几何感知**（不使用 `edge_lengths`/`face_areas`/level set 分数）。

### 1.3 EB（嵌入式边界）

| 路径 | 语义 | 证据 |
|---|---|---|
| 非 ECT（Yee/CKC/PSATD/hybrid/macroscopic） | **阶梯掩码**：被覆盖或相邻的格点不更新 | `EmbeddedBoundaryInit.cpp`（`MarkUpdateCellsStairCase`），初始化为 `Source/Initialization/WarpXInitData.cpp:1640+` |
| ECT | **cut-cell/conformal**：`edge_lengths`、`face_areas`、`area_mod` + borrowing；语义是 **PEC**（覆盖侧切向 E=0） | 初始化 `WarpXInitData.cpp:1607-1632`；B 更新 `Source/FieldSolver/FiniteDifferenceSolver/EvolveB.cpp:117,228`（`EvolveBCartesianECT`，注释限定 2D3V/3D3V）；面积扩展 `WarpXFaceExtensions.cpp` |
| EB 体是否为介质 | **否**：`eb_update_E=0` 把 EB 内的场冻结；没有“介质内部继续求解”的路径 | 同上；旧调研文档 §0.2/§3.4 亦同结论 |

**结论：当前没有“介质 + cut-cell”能力**。`ECT` 是 PEC 语义；`macroscopic` 是阶梯；二者互斥。
这正是 C4（介质球/柱）、C6（亚网格薄层）、C10（核壳球）的阻塞项（记为 M2）。

### 1.4 源与诊断

| 能力 | 实测 | 备注 |
|---|---|---|
| 激光/平面波注入 | `lasers.names` + `<laser>.position/direction/polarization/e_max/wavelength/profile`；`profile = parse_field_function(X,Y,t)` | 见 `Docs/source/usage/parameters.rst:2435-2540`；**任意 $k_x$ 相位斜坡可用 `field_function` 实现**，但需 G0 验证 |
| TF/SF 面源 | **无**（上游与本 fork 都没有 `PlaneWaveSource`） | 旧树里的 `PlaneWaveSource` 是自研代码，已随树删除 |
| 场输出 | plotfile（AMReX） | openPMD 当前 OFF |
| 点/线/面探针 | `FieldProbe`（`Point/Line/Plane`） | `Source/Diagnostics/ReducedDiags/FieldProbe.H:25-31` |
| 能量/坡印廷/表面电荷 | `FieldEnergy`、`FieldPoyntingFlux`、`ChargeOnEB` | `Source/Diagnostics/ReducedDiags/MultiReducedDiags.cpp:59-73` |
| 近远场变换（NTFF） | **无** | 散射截面用**近场闭合面坡印廷通量**替代 |
| PML | `boundary.field_lo/hi = pml` | `parameters.rst:1018-1032` |

### 1.5 流体/矩求解器（本项目新增部分）

- `Source/Fluids/`：`WarpXFluidContainer`（冷相对论流体：密度 + 动量，MUSCL-Hancock + SSP-RK3）。
- `Source/MomentsSolver/`（**2026-10-31 新增，尚未接入 CMake**）：
  - 方程对象：`Euler5`（5 矩 Euler，$\gamma$ 可配）、`TenMoment`（10 矩压力张量，`k0` 闭合系数）、
    `Maxwell6`（真空 6 分量）；
  - 1D 有限体积驱动 `FV1D`：MC 限制器 + LLF 数值通量 + SSP-RK3；周期/出流边界；
  - 已有单元测试（`Tests/`）：Sod 激波管守恒、10 矩各向异性平衡保持、Maxwell 脉冲以 $c$ 传播；
  - **没有**：多维、磁场洛伦兹源项、与 Maxwell 的耦合、与 WarpX 场求解器的接口。
- 结论：矩求解器目前是"1D 方程骨架 + 单元测试"，本计划的 C8/C9 用来标定它将来要接的
  **等离子体色散与磁化波模**（先用 PIC/混合求解器给出基线）。

### 1.6 能力矩阵（C1–C10 × 现在能不能跑）

图例：✅ 现在可跑；🟡 需小改（构建/脚本）；🔴 需新功能（阻塞）。

| 算例 | 现在 | 阻塞项 |
|---|---|---|
| C1 | ✅（需 1D/2D 构建） | — |
| C2-0（PEC 极限自检） | ✅（ECT；需 2D/3D 构建） | — |
| C2（介质 cut-cell） | 🔴 | M2：介质 + EB 场更新 |
| C3 | ✅ | — |
| C4-PEC | ✅（ECT 球/柱） | — |
| C4-介质 | 🔴 | M2：介质 + EB 的场更新 |
| C5 | ✅ | — |
| C6-厚层 | ✅ | — |
| C6-薄层（$d<\Delta x$） | 🔴 | M2：亚网格界面平均 |
| C7 | ✅ | — |
| C8 | ✅（PIC / 冷流体） | — |
| C9 | 🟡（示例是 PICMI，需 PYTHON=ON 或改写为 inputs） | — |
| C10a | 🔴 | M2 + M4（AMR 材料） |
| C10b | 🔴 | M2 + M4 + 磁化等离子体耦合 |

---

## 2. 方法论：四道门

### 2.1 门 G0：文献卡（跑之前）

每个算例提交 `LITERATURE.md`，必须包含：

1. 公式出处：作者、期刊、年、卷、页、式号、DOI（本文件 §9 已给出可用条目）；
2. 约定：时间因子、$\varepsilon$ 虚部符号、$\mu_r$、归一化（$c=1$ 还是 SI）；
3. **要复现的数值表**：例如 Mie 的 $Q_{\rm ext}(x)$ 表、Airy 曲线的极值位置、MacLeod 书中的
   滤光片设计例；
4. 适用条件与失效边界（例如：Fresnel 公式要求界面无限大、平面、均匀；Mie 要求球、均匀、各向同性）。

**G0 的通过判据**：`theory.py` 能在 3 个独立算例上复现文献数值到文献给的有效位数（通常 3–4 位）。

### 2.2 门 G1：解析解独立实现（两套）

- 一套**手写**（numpy 直接实现公式，禁止调用现成库的“黑盒”）；
- 一套**独立来源**（第三方库或文献表格），逐点对比；
- 两套不一致时，以文献表格为准，并记录差异原因；
- 所有 `theory.py` 必须有单元测试：`tools/test_theory.py`（pytest），随仓库提交。

### 2.3 门 G2：数值实验设计

每个算例固定做四件事（缺一不可）：

1. **真空对照（仪器地板）**：删掉界面/物体（或令 $\varepsilon_r=1$、$\sigma=0$），其余参数不变，
   测同一套分析量。它给出该算例的"数值噪声本底"。**所有物理判据必须显著高于本底**。
2. **分辨率扫描**：至少 3 个 $\Delta x$（例如 $\lambda/20,\lambda/40,\lambda/80$），报告收敛阶。
3. **PML 距离扫描**：PML 距界面/物体至少 2 个波长起，验证结果对 PML 位置不敏感（<0.2%）。
4. **时间窗**：记录量在"波前已通过、尚未到域边界"的窗口内取值；不用全局能量判稳定性。

### 2.4 门 G3：数据完整性（防造假的具体机制）

| 机制 | 做法 |
|---|---|
| 原始数据只读 | `raw/` 目录权限 `chmod -R a-w`，只在 `run.sh` 里写入一次；分析只读不写 |
| 溯源清单 | `meta.json`：二进制 sha256、输入 sha256、`CMakeCache.txt` 关键项、编译器版本、主机名、时间、命令 |
| 单一数据源 | 报告与图表只从 `results.json` 生成；`results.json` 只由 `analyze.py` 从 `raw/` 生成 |
| 重算校验 | `tools/verify.py`：重跑 `analyze.py`，逐字段比对 `results.json`，不一致即失败 |
| 预注册 | `expected.json`（阈值 + 判据）必须在**第一次运行前**提交，之后只能通过 `CHANGELOG.md` 修改并说明理由 |
| 失败留档 | 失败算例保留 `raw/` 与 `FAILED.md`；`STATUS.md` 汇总 |
| 版本 | `WarpX_fluid` 当前 git **无任何提交**：溯源必须用**文件 sha256 清单**（`Source/` 全量 `sha256sum`）而不是 commit hash |
| 禁改数据 | 明文规定：任何对 `raw/`、`results.json` 的手工编辑都视为数据造假，算例作废并全量重跑 |

---

## 3. 通用数值约定

- **单位**：SI，$\lambda_0=1\ \mathrm{m}$，$f_0=c/\lambda_0=299.792458\ \mathrm{MHz}$；时间以周期
  $T_0=1/f_0$ 计。所有长度以 $\lambda_0$ 归一化后报告。
- **CFL**：Yee 取 0.5（3D 稳定上限 $1/\sqrt3\approx0.577$）；ECT 取 0.4（面积扩展后局部 CFL 更严）。
- **PML**：12–16 格，$\sigma$ 多项式 3 阶；PML 内不放置任何界面/物体。
- **分辨率**：平面界面 $\ge 20$ 格/$\lambda_{\rm med}$；cut-cell 曲面 $\ge 30$ 格/$\lambda$（且 $\ge a/20$，
  曲率半径 $\ge 10$ 格）；薄层问题按 §C6 单独规定。
- **源**：脉冲优先（宽带，一次运行覆盖一个频段），Gaussian 时间包络，中心频率 $f_0$，
  带宽 $\ge 60\%$。用激光注入时记录 `laser.e_max` 与到达界面的幅度。
- **采样**：每个场的监测面/线/点以 `FieldProbe` 输出；时间步采样 $\le T_0/40$。
- **窗口**：所有拟合窗必须比波前位置至少提前 2 个 $\lambda$（避免未填充域伪影）。

---

## 4. 十个算例

> 每个算例的模板：目标 → 依赖 → 物理与理论解 → 数值设置 → 测量方法 → 判据 → 交付物 →
> 风险与对策 → 文献门槛。

---

### C1　平面介质界面：Fresnel 反射/折射与 Brewster 角

**目标**：在**最干净**的几何（无限平面、均匀介质）上标定材料求解器的 $r$、$t$、$\theta_t$、能量
守恒，并测出 p 偏振 Brewster 角。

**依赖**：无（1D 正入射段不需要 EB；斜入射段用"倾斜界面"实现，见下）。需 1D 与 2D 构建。

**物理与理论解**（完整推导见附录 A.1）
- 入射：$\mathbf{k}_i=k_0 n_1(\sin\theta_i,0,\cos\theta_i)$；Snell：$k_x$ 连续，
  $k_{zj}=k_0\sqrt{n_j^2-n_1^2\sin^2\theta_i}$（取 $\mathrm{Re}>0,\mathrm{Im}\ge0$ 分支）。
- s 偏振：$r_s=\dfrac{\mu_2k_{z1}-\mu_1k_{z2}}{\mu_2k_{z1}+\mu_1k_{z2}}$，
  $t_s=\dfrac{2\mu_2k_{z1}}{\mu_2k_{z1}+\mu_1k_{z2}}$。
- p 偏振：$r_p=\dfrac{\varepsilon_2k_{z1}-\varepsilon_1k_{z2}}{\varepsilon_2k_{z1}+\varepsilon_1k_{z2}}$，
  $t_p=\dfrac{2n_1\cos\theta_i}{n_2\cos\theta_i+n_1\cos\theta_t}$（E 振幅约定，$\mu_r=1$）。
- 功率：$R=|r|^2$，$T=\dfrac{\mathrm{Re}\,k_{z2}}{\mathrm{Re}\,k_{z1}}|t|^2$，无耗时 $R+T=1$。
- Brewster：$\tan\theta_B=n_2/n_1$（$\mu_r=1$、无耗），此时 $r_p=0$。

**数值设置**
- *子算例 1a（1D 正入射）*：`WARPX_DIM_1D_Z`，$\varepsilon_r\in\{2.25,4,9,16\}$，$\sigma=0$，
  界面在 $z=0$，真空侧 $z<0$；PML 两侧；脉冲注入；两平面采样做 $+z/-z$ 分解。
- *子算例 1b（2D 斜入射）*：`WARPX_DIM_XZ`，界面相对 $z$ 轴倾斜 $\alpha\in\{0,15,30,45,60,75\}^\circ$；
  入射沿 $z$（法向源），界面倾斜 ⇒ 等价于入射角 $\theta_i=\alpha$ 打平面界面；
  侧向 PML；反射波沿 $2\alpha$ 方向传播，测量窗取在反射方向上。
- *子算例 1c（有耗）*：$\varepsilon_r=4$，$\sigma\in\{0.01,0.1,1\}$ S/m，$f=1$ GHz，复 Fresnel。
- *子算例 1d（p 偏振/Brewster）*：$\varepsilon_r=4$（$\theta_B=63.43^\circ$）、$\varepsilon_r=9$
  （$\theta_B=71.57^\circ$），在 $\theta_B\pm5^\circ$ 内扫 5 个角度。
- *子算例 1e（仪器地板）*：同几何、$\varepsilon_r=1$（无界面），测同一分析量的本底。

**测量方法**
- 场分解：在真空区取两条 $z=$ 常数平面，时间锁相拟合
  $E_y(z,t)=a\cos(\omega t-kz)+b\cos(\omega t+kz)+c$ ⇒ $A=a-ib$（复振幅），$r=A_{\rm ref}/A_{\rm inc}$。
  （这是对旧报告 P1 方法的复用；解析解与拟合都在 `analyze.py` 里实现并单测。）
- 透射：介质侧同样做双波拟合，得 $t$ 与 $k_{\rm med}$；$\theta_t=\arcsin(k_x/k_{\rm med})$。
- 能量：$R+T$ 与 $1$ 的偏差；所有量在波前通过后的窗口内取。

**判据（预注册）**

| 量 | A 级（必须） | B 级（目标） |
|---|---|---|
| $\lvert r\rvert,\lvert t\rvert$ | $\le 2\%$ | $\le 1\%$ |
| $\arg r$（相位） | $\le 4^\circ$ | $\le 2^\circ$ |
| $\theta_t$ | $\le 0.8^\circ$ | $\le 0.5^\circ$ |
| $R+T$（无耗） | $1\pm1\%$ | $1\pm0.5\%$ |
| Brewster：$\lvert r_p\rvert$ 在 $\theta_B$ | $\le 0.02$ | $\le 0.01$ |
| 仪器地板（1e） | $<10^{-4}$（相对入射） | $<10^{-5}$ |

**交付物**：`cases/c01_fresnel/`（inputs、gen.py、theory.py、analyze.py、expected.json、run.sh、
meta.json、raw/、results.json、figs/、LITERATURE.md、README.md、sha256sums.txt）。

**风险与对策**
- 1D/2D 构建未验证材料路径（代码注释说只支持 3D）→ 1a 先跑 1D，若失败则退到 3D 平板（薄板）。
- 倾斜界面与 PML 相交处的边角散射 → 侧向 PML 距界面端点 ≥2λ；用 1e 的本底扣除。
- 激光源的谱纯度 → 用宽带脉冲 + FFT，不做单频驱动。

**文献门槛（G0）**：Born & Wolf《Principles of Optics》§1.5（Fresnel 公式）；Jackson
《Classical Electrodynamics》§7.3；有耗介质见 Bohren & Huffman §2（复折射率约定）。
本算例的 `theory.py` 必须复现 Born & Wolf §1.5 的 $R_s,R_p$ 曲线（3 位有效数字）。

---

### C2　斜平面介质界面 cut-cell：C1 的「几何实现对照」

**定位（与 C1 的关系）**：C2 **不是新物理**，而是 **C1 的逐参数复制**：源、域、PML、监测面、
参数扫描、测量方法、判据全部沿用 C1；**唯一差别**是界面的表示方式——

- **C1**：界面由 `macroscopic.epsilon_function(x,z)` 以**阶梯**方式给出（网格内 ε 分片常数）；
- **C2**：界面由 `warpx.eb_implicit_function` 定义，由 **cut-cell 几何**（level set +
  `edge_lengths`/`face_areas`/`area_mod`）表示，材料系数由界面分数做**亚网格平均**（subpixel）。

因此 **C2 的净效应 = C2 结果 − C1 结果**：把「物理误差（Fresnel/Snell/损耗）」与
「几何误差（界面位置、面元、平均）」分开。这是本算例存在的唯一理由。

**依赖**：M2（介质 + EB 的场更新）。M2 完成前只能跑 **C2-0**。

**C2-0（现在即可跑）**：把同一倾斜界面换成 **PEC**（`algo.maxwell_solver = ect`），
理论解 $r=-1$。它不验证介质，但给出 cut-cell **几何误差的下界**（相位、伪散射、能量），
并作为 C2 的「几何本底」。C2-0 的判据在表中单列。

**物理与理论解**（与 C1 完全一致，见附录 A.1）
- 入射沿 $z$、界面法向 $\mathbf n=(\sin\alpha,0,\cos\alpha)$ ⇒ 入射角 $\theta_i=\alpha$；
- 复 Fresnel $r_s,r_p,t_s,t_p$、Snell、Brewster $\tan\theta_B=n_2/n_1$；
- 反射波方向 $2\alpha$，透射波方向由 Snell 与 $\mathbf n$ 决定；
- **额外要求**：亚网格界面的「有效位置」应等于解析平面（偏差 $<\Delta/10$ 才算实现正确），
  由反射相位反推界面位置来检验。

**参数扫描**（角度从小到大、ε 变化、偏振、损耗、亚网格位置）

| 扫描 | 取值 | 子算例 | 目的 |
|---|---|---|---|
| 入射角 $\alpha$ | 0°, 5°, 10°, 15°, 30°, 45°, 60°, 75°, 80° | 2a | 小角度≈网格对齐、大角度≈掠入射；两端最易暴露问题 |
| $\varepsilon_r$（无耗） | 1.5, 2.25, 4, 9, 16 | 2b | 弱→强折射率；Brewster 角随之变化 |
| $\sigma$（$\varepsilon_r=4$） | 0.01, 0.1, 1 S/m | 2c | 损耗 + cut-cell 的耦合 |
| 偏振 | s、p | 2a/2d | Brewster 是否仍消失 |
| 分辨率 | $\lambda/20,\lambda/40,\lambda/80$ | 2e | 收敛阶（阶梯一阶；subpixel 期望接近二阶） |
| 亚网格偏移 | 界面相对网格面偏移 0, ¼, ½, ¾ $\Delta$ | 2f | cut-cell 对界面位置的敏感性（C1 做不到） |
| 几何极限 | $\alpha\to0$ 必须回到 C1a；$\alpha=80^\circ$ 掠入射 | 2g | 一致性与极端角 |
| 真空对照 | $\varepsilon_r=1$ | 2h | 仪器地板 |
| PEC 极限 | 同几何 + ECT | 2i | C2-0 几何本底 |

**数值设置**
- 2D `XZ`，域 $40\lambda\times40\lambda$，PML 16 格；CFL=0.4（ECT）/0.5（材料）；
- 主分辨率 $\lambda/40$；2e 用 $\lambda/20,\lambda/40,\lambda/80$；
- 时间 30 个周期；所有量在「波前已通过、未到 PML」的窗口内取；
- **每个 C2 子算例必须与 C1 的同参数运行成对出现**（同域、同源、同监测面）。

**测量方法**
- 与 C1 相同：双波分解、双波拟合、能量；
- **几何误差**：$E_{\rm geo}=C2-C1$（同参数逐点），报告 L2 与相位差；
- **解析误差**：C2 对 Fresnel 的偏差（判据同 C1）；
- **界面位置**：由反射相位反推有效界面位置，与解析平面比较；
- **能量**：介质 + 真空总场能（不含 PML）在反射/透射完成后守恒。

**判据（预注册）**

| 量 | C2（介质 cut-cell，A 级） | C2-0（PEC 极限，现在可跑） |
|---|---|---|
| 对 Fresnel 的 $\lvert r\rvert,\lvert t\rvert$ | $\le 2\%$ | $\lvert r+1\rvert\le 1\%$ |
| 相位 | $\le 4^\circ$ | $\le 2^\circ$ |
| $\theta_t$ | $\le 0.8^\circ$ | — |
| 能量 | $1\pm1\%$ | $1\pm10^{-3}$ |
| 与 C1 同参数之差 | $\le1\%$（幅度）、$\le2^\circ$（相位） | — |
| 亚网格偏移 $0\to\tfrac12\Delta$ | $\le1\%$ | $\le1\%$ |
| 收敛阶 | $\ge1.5$（subpixel 期望接近 2） | $\ge1.5$ |
| 掠入射 $75^\circ$–$80^\circ$ | $\le5\%$（单列，不进主判据） | 单列 |
| 仪器地板 | $<10^{-4}$ | 同 |

**交付物**：`cases/c02_cutcell_dielectric/`：`inputs`、`gen.py`、`theory.py`、`analyze.py`、
`expected.json`、**与 C1 的配对表 `delta_vs_c1.csv`**、收敛图、偏移敏感性图、C2-0 子目录、
`LITERATURE.md`、`README.md`、`sha256sums.txt`。

**风险与对策**
- 掠入射时反射波几乎与界面平行，测量窗会与 PML 重叠 → 单列、不作主判据；
- 若 M2 用「EB 内常数 ε」（阶跃）而非 subpixel，预期一阶收敛；本算例正是量化这个差别——
  它决定 M2 是否必须做亚网格平均；
- C1 自身的阶梯误差会污染 $\Delta$：C1 必须先达到 B 级；
- `eb_implicit_function` 的符号决定「哪一侧是内部」，必须与材料侧一致，否则界面两侧 ε 会反
  （旧报告曾在这条上踩过坑）；
- ECT 与材料 Yee 的语义不同：C2 的介质场更新必须走「材料 Yee」，不能把 ECT 的 PEC 公式
  直接改成介质版（会破坏两种语义）。

**文献门槛（G0）**：C1 的文献（Born & Wolf §1.5；Jackson 第 7 章）+
Dey & Mittra 1997（cut-cell 几何，DOI 10.1109/75.622536）+
Farjadpour 2006（subpixel，DOI 10.1364/OL.31.002972）+
Deinega & Valuev 2007（有耗 subpixel，DOI 10.1364/OL.32.003429）+
Xiao & Liu 2004（enlarged cells，DOI 10.1109/LMWC.2004.837384，WarpX `area_mod` 即此思想）。
**要复现**：Farjadpour 文中斜平面/曲面界面的反射误差随分辨率的曲线（以论文为准，
G1 复核后再把数值写死）。

---

### C3　介质平板：Fabry–Perot 干涉

**目标**：验证多次反射、薄层相位、以及"厚度 < 网格"时的误差行为；给出透射谱与 Airy 公式的
逐频点对比。

**依赖**：材料求解器（现成）。斜入射段用倾斜平板（同 C1b 的几何技巧）。

**物理与理论解**（附录 A.2）
- 单层板（$n_2$，厚 $d$，两侧 $n_1$）：$r_{12}=\dfrac{n_1-n_2}{n_1+n_2}$，$\delta=k_0n_2d\cos\theta_t$，
  $$r=\frac{r_{12}+r_{23}e^{2i\delta}}{1+r_{12}r_{23}e^{2i\delta}},\qquad
    t=\frac{t_{12}t_{23}e^{i\delta}}{1+r_{12}r_{23}e^{2i\delta}}$$
- 无耗时 $T=\dfrac{(1-R)^2}{(1-R)^2+4R\sin^2\delta}$，$R=|r_{12}|^2$；
  自由光谱范围 $\mathrm{FSR}=\dfrac{c}{2n_2d\cos\theta_t}$，精细度 $\mathcal F=\dfrac{\pi\sqrt R}{1-R}$。
- 极值：$\delta=m\pi$ 时 $T_{\max}=1$（$R$ 任意）；$\delta=(m+\tfrac12)\pi$ 时
  $T_{\min}=\left(\dfrac{1-R}{1+R}\right)^2$。

**数值设置**
- 3a（1D 正入射）：$n_2\in\{1.5,2,3\}$，$d\in\{\lambda_0/2,\lambda_0,2\lambda_0\}$；
  宽带脉冲 0.5–1.5 $f_0$；测 $T(f)$。
- 3b（2D 斜入射）：平板倾斜 $\alpha=30^\circ,45^\circ$；测 $T(f)$ 与 $\mathrm{FSR}$ 的角度依赖
  （$\cos\theta_t$ 因子）。
- 3c（亚网格）：$d=\lambda_0/20$（$<\Delta x$ 的一半），同样的 $T(f)$；**预期失败**，
  用于量化 M2（亚网格平均）的必要性。
- 3d（仪器地板）：去掉平板。

**测量方法**
- $T(f)=|\tilde E_{\rm trans}(f)|^2/|\tilde E_{\rm inc}(f)|^2$，FFT 加 Hann 窗；
  入射谱用 1e 的同一次运行或同一源的参考运行。
- 峰位用抛物线拟合；FSR = 相邻峰间隔；finesse = FSR/半高全宽。

**判据（预注册）**

| 量 | 厚层（$d\ge\lambda/2$） | 亚网格（3c） |
|---|---|---|
| $T(f)$ 逐点 | $\le 2\%$ | 记录（预期 >10%，作为 M2 的动机） |
| 峰位 | $\le 0.5\%$ | — |
| FSR | $\le 1\%$ | — |
| finesse | $\le 5\%$ | — |
| $T_{\max}$ | $\ge 0.98$（无耗） | — |
| 能量 $R+T$ | $1\pm0.5\%$ | — |

**交付物**：`cases/c03_fabry_perot/`（含 $T(f)$ 曲线、峰位表、亚网格失败记录）。
**风险**：脉冲带宽不足导致峰分辨不清 → 带宽 ≥60%，并用 3 个带宽做敏感性检查。
**文献门槛**：MacLeod《Thin-Film Optical Filters》第 2 章（Airy 公式与 FSR/finesse 定义，
DOI 10.1201/9781420073034）；Born & Wolf §7.6（多光束干涉）。**要复现**：MacLeod 书中
单层滤光片的 $T(\lambda)$ 图（3 位）。

---

### C4　圆柱与球的 Mie 散射

**目标**：把验证从"平面"推进到"曲面"，覆盖 2D（无限圆柱）与 3D（球）、PEC 与介质两种边界；
输出散射截面 $Q_{\rm ext}$、角分布 $S(\theta)$、近场误差。

**依赖**：PEC 段用 ECT（现成）；介质段需 M2。

**物理与理论解**（附录 A.4/A.5）
- 球（半径 $a$，$x=ka$，$m=n_s/n_m$）：
  $$a_n=\frac{m\psi_n(mx)\psi_n'(x)-\psi_n(x)\psi_n'(mx)}{m\psi_n(mx)\xi_n'(x)-\xi_n(x)\psi_n'(mx)},\quad
    b_n=\frac{\psi_n(mx)\psi_n'(x)-m\psi_n(x)\psi_n'(mx)}{\psi_n(mx)\xi_n'(x)-m\xi_n(x)\psi_n'(mx)}$$
  $$Q_{\rm ext}=\frac{2}{x^2}\sum_{n=1}^{\infty}(2n+1)\,\mathrm{Re}(a_n+b_n),\quad
    Q_{\rm sca}=\frac{2}{x^2}\sum(2n+1)(|a_n|^2+|b_n|^2)$$
  PEC 极限：$a_n=\psi_n'(x)/\xi_n'(x)$，$b_n=\psi_n(x)/\xi_n(x)$。
- 圆柱（半径 $a$，$x=ka$，正入射）：TM（$E\parallel$ 轴）与 TE 各一组系数（附录 A.5），
  $Q_{\rm ext}^{\rm TM}=\frac{2}{x}\sum_{n=-\infty}^{\infty}\mathrm{Re}(a_n)$（2D 归一化约定在 G1 与
  Bohren & Huffman 表核对后固化）。
- 近场：$\mathbf E_{\rm sca}$ 的级数解在 $r=1.5a$ 的圆/球面上逐点可算，用于**无需 NTFF** 的对比。

**数值设置**
- 4a（2D 圆柱，PEC）：$x\in\{0.5,1,2,5,10\}$；$a=x\lambda_0/(2\pi)$；$\Delta x\le\min(a/20,\ \lambda_0/30)$；TM 与 TE；PML 距柱心 ≥5λ。
- 4b（2D 圆柱，介质）：$m=1.5$，同上；需 M2（先做阶梯基线，记录误差）。
- 4c（3D 球，PEC）：$x\in\{0.5,1,2,5\}$；$\Delta x\le\min(a/20,\ \lambda_0/30)$；二阶收敛检查。
- 4d（3D 球，介质）：$m=1.5,2$；M2 后做。
- 4e（仪器地板）：无物体。

**测量方法**
- $Q_{\rm ext}$：在包围物体的闭合面（2D 为圆、3D 为球）上积分
  $\oint \mathrm{Re}(\mathbf E\times\mathbf H^*)\cdot d\mathbf S$ / 入射强度 / 几何截面。
- $S(\theta)$：远场用近场的 $\mathbf E_{\rm sca}$ 在 $r=1.5a$ 上做傅里叶–贝塞尔展开
  （2D）或矢量球谐投影（3D），避免实现 NTFF。
- 近场误差：$\|\mathbf E_{\rm sim}-\mathbf E_{\rm Mie}\|_2/\|\mathbf E_{\rm Mie}\|_2$，
  在物体外 0.5λ 的壳层内。

**判据（预注册）**

| 量 | 2D | 3D |
|---|---|---|
| $Q_{\rm ext}$（$x\ge1$） | $\le 3\%$ | $\le 2\%$ |
| $Q_{\rm ext}$（$x<1$） | $\le 5\%$ | $\le 5\%$ |
| 近场 L2 | $\le 5\%$ | $\le 3\%$ |
| $S(\theta)$ 角分布 rms | $\le 8\%$ | $\le 5\%$ |
| 收敛阶 | $\ge1.5$ | $\ge1.5$ |
| 仪器地板 | $Q_{\rm ext}<10^{-3}$ | 同 |

**交付物**：`cases/c04_mie/`（$Q_{\rm ext}(x)$ 表、角分布、近场切片、收敛图、
PEC/介质/阶梯三套对照）。
**风险**：
- 阶梯球在低分辨率下 $Q_{\rm ext}$ 系统性偏高 → 用 4 个分辨率外推，并单列"阶梯偏差"曲线；
- 2D 圆柱的 $Q_{\rm ext}$ 归一化常数易错 → G1 用文献表逐个 $x$ 核对后才允许进结果；
- 曲面在 ECT 下的最小面分数过小 → 记录 `area_mod` 统计，必要时加局部加密（C10）。

**文献门槛**：Mie 1908（DOI 10.1002/andp.19083300302）；Bohren & Huffman
（DOI 10.1002/9783527618156；球 Mie 级数与无限柱级数见该书对应章节）；Wiscombe 1980
（DOI 10.1364/AO.19.001505，算法与表格）；van de Hulst《Light Scattering by Small Particles》
（柱/球级数的经典写法）。**要复现**：Bohren & Huffman 书中 $m=1.5$ 的 $Q_{\rm ext}(x)$ 表
（$x=0.1$–$10$，4 位有效数字）。

---

### C5　有耗介质：复折射率、趋肤深度、吸收率

**目标**：验证 $\sigma$ 进入 $\alpha/\beta$ 系数的方式、损耗介质的 $R/T$、以及趋肤深度。

**依赖**：材料求解器（现成，`macroscopic.sigma`）。

**物理与理论解**（附录 A.1/A.6）
- 复介电常数 $\tilde\varepsilon_r=\varepsilon_r-i\dfrac{\sigma}{\omega\varepsilon_0}$，
  复折射率 $\tilde n=\sqrt{\tilde\varepsilon_r\mu_r}$。
- 平面波 $e^{i(kz-\omega t)}$：$k=k_0\tilde n$，场 $\propto e^{-\kappa z}$，
  $\kappa=\mathrm{Im}k$，$\delta_s=1/\kappa$。
- 良导体极限：$\delta_s=\sqrt{2/(\omega\mu\sigma)}$。
- 界面（正入射）：$r=\dfrac{\tilde n_1-\tilde n_2}{\tilde n_1+\tilde n_2}$；
  吸收率 $A=1-R-T$。

**数值设置**
- 1D：$\varepsilon_r=4$，$\sigma\in\{0.01,0.1,1,10\}$ S/m，$f=1$ GHz，
  $\delta_s$ 从 $0.53\lambda$ 到 $0.017\lambda$（跨约 30 倍）；分辨率需 $\Delta x\le\delta_s/8$
  （对 $\sigma=10$ 需 $\Delta x\le\lambda/480$ —— **这是本算例的刻意压力点**，用来暴露
  "薄趋肤层 + 均匀网格"的代价，并作为 AMR 动机记录）。
- 2D：倾斜 30° 的有耗平板（同 C1b 几何）。
- 仪器地板：$\sigma=0$ 同几何。

**测量方法**
- 透射场在介质内两个深度 $z_1,z_2$ 的幅度比 ⇒ $\kappa$；
- $R$ 用真空侧双波拟合；$T$ 用介质后侧单波拟合（考虑衰减修正）；
- $A=1-R-T$；与 $\delta_s$ 公式对比。

**判据**：$\kappa$（或 $\delta_s$）$\le 3\%$；$R,T\le 2\%$；$R+T+A=1\pm1\%$；
良导体极限（$\sigma=10$）$\delta_s\le 5\%$。

**交付物**：`cases/c05_lossy/`（$\delta_s$ 表、$R,T,A$ 表、压力点记录）。
**风险**：$\sigma=10$ S/m 时 $\delta_s\ll\Delta x$，误差不可避免 → 明确记录"该点需要 AMR/亚网格"，
并把它作为 C10/M4 的输入。
**文献门槛**：Jackson《Classical Electrodynamics》第 7–8 章（反射/折射；导体中的衰减与趋肤深度）；Bohren & Huffman §2（复折射率约定与符号）。

---

### C6　多层介质：增透膜与 Bragg 反射镜

**目标**：验证多界面叠加、相位累积、以及"层厚 $<\Delta x$"的亚网格能力；给出特征矩阵逐频点对比。

**依赖**：材料求解器（厚层现成）；薄层需 M2。

**物理与理论解**（附录 A.3）
- 特征矩阵：$M_j=\begin{pmatrix}\cos\delta_j & \frac{i}{\eta_j}\sin\delta_j\\
  i\eta_j\sin\delta_j & \cos\delta_j\end{pmatrix}$，$\delta_j=k_0n_jd_j\cos\theta_j$，
  $\eta_j=n_j\cos\theta_j$（s）或 $n_j/\cos\theta_j$（p）（$\mu_r=1$）。
- $[B;C]=\left(\prod M_j\right)\begin{pmatrix}1\\ \eta_s\end{pmatrix}$，
  $r=\dfrac{\eta_0B-C}{\eta_0B+C}$，$t=\dfrac{2\eta_0}{\eta_0B+C}$，$T=\dfrac{\eta_s}{\eta_0}|t|^2$。
- 单层增透：$n_c=\sqrt{n_0n_s}$，$d=\lambda_0/(4n_c)$ ⇒ $R(\lambda_0)=0$。
- 四分之一波长 Bragg 堆（$N$ 对，$n_0=n_s=1$）：
  $R_{\max}=\left(\dfrac{1-q}{1+q}\right)^2$，$q=(n_H/n_L)^{2N}$。

**数值设置**
- 6a（1D 增透膜）：$n_s=1.5$，$n_c=1.2247$，$d=\lambda_0/(4n_c)$；宽带 0.7–1.3 $f_0$；测 $R(f)$。
- 6b（1D Bragg）：$n_H=2.3$，$n_L=1.45$，$N=10$；测阻带中心与宽度、$R_{\max}$。
- 6c（斜入射）：6a/6b 倾斜 30°（s/p 分别）。
- 6d（亚网格）：把 6a 的膜厚改为 $d=\lambda_0/20$；预期失败，量化 M2 需求。
- 仪器地板：无膜。

**测量方法**：$R(f)$ 由真空侧反射脉冲的 FFT 得到（须扣除源谱）；阻带中心 = $R>0.9$ 区间中点；
与特征矩阵逐频点比较。

**判据**：厚层 $R(f)$ 逐点 $\le2\%$；$R_{\min}$（增透）$\le0.5\%$；Bragg $R_{\max}\le2\%$、
中心频率 $\le1\%$；亚网格记录（不设通过线，作为 M2 证据）。

**交付物**：`cases/c06_multilayer/`（设计参数、$R(f)$、阻带表、亚网格记录）。
**风险**：宽带脉冲在薄膜内多次反射导致时域重叠 → 用足够长的窗口与 $\ge$ 20 个周期。
**文献门槛**：MacLeod 第 2–3 章（特征矩阵、AR 设计、Bragg 堆）；
Born & Wolf §7.6。**要复现**：MacLeod 中 10 层 $\mathrm{TiO_2/SiO_2}$ 滤光片的 $R(\lambda)$ 曲线。

---

### C7　全内反射与倏逝波

**目标**：验证倏逝场的指数衰减、$R=1$、以及 PML 对倏逝场的吸收（不引入伪反射）。

**依赖**：材料 + PML（现成）。

**物理与理论解**（附录 A.6）
- $\sin\theta_c=n_2/n_1$（$n_1>n_2$）。
- $\theta_i>\theta_c$ 时 $k_{z2}=i\kappa$，$\kappa=k_0\sqrt{n_1^2\sin^2\theta_i-n_2^2}$，
  介质外场 $\propto e^{-\kappa z}$，衰减长度 $1/\kappa$。
- $R=1$（无耗）；相位 $\phi$ 随角度变化；Goos–Hänchen 位移
  $D=-\dfrac{1}{k_0n_1\cos\theta_i}\dfrac{d\phi}{d\theta_i}$（选做，实施前在 G0 复核符号与约定）。

**数值设置**
- 2D：$n_1=1.5$（介质），$n_2=1$（真空），$\theta_c=41.81^\circ$；
  $\alpha\in\{45,50,60,70\}^\circ$；PML 在真空侧；探针沿法向取 20 点。
- 仪器地板：$\alpha=30^\circ<\theta_c$（应有透射，作为对照）。

**测量方法**：$R$ 用双波分解；$\kappa$ 用 $\ln|E|$ 对法向距离的线性拟合（$r^2>0.999$）；
PML 反射用"远离界面处的反向波"幅度估计。

**判据**：$R=1\pm0.1\%$；$\kappa$（或 $1/\kappa$）$\le3\%$；PML 伪反射 $<10^{-3}$。
**交付物**：`cases/c07_tir/`（衰减曲线、$\kappa(\theta)$ 表）。
**风险**：倏逝场在 PML 内被放大（PML 对倏逝波不是严格吸收）→ 记录 PML 厚度扫描。

**文献门槛**：Born & Wolf §1.5（全反射与倏逝波）；Jackson 第 7 章；GH 位移见 Artmann 1948
（在 G0 中补 DOI）。

---

### C8　无磁化等离子体：色散、截止、趋肤

**目标**：把验证从"无源介质"推进到"等离子体"：这是磁流体代码将来要耦合的对象。
同时是"粒子路径（PIC）"与"流体路径（冷流体）"的交叉验证。

**依赖**：PIC（现成）与 `WarpXFluidContainer`（现成）。

**物理与理论解**（附录 A.7）
- $\omega_p=\sqrt{n_ee^2/(\varepsilon_0m_e)}$；
  $\omega^2=\omega_p^2+c^2k^2$（无磁化、冷等离子体）；
- 截止：$\omega<\omega_p$ 时 $k$ 为虚，$\kappa=\sqrt{\omega_p^2-\omega^2}/c$，趋肤深度 $1/\kappa$；
- 半空间反射（正入射）：$r=\dfrac{k_0-k}{k_0+k}$（复 $k$）。

**数值设置**
- 1D：$n_e=10^{18}\ \mathrm{m^{-3}}$（$f_p\approx8.98$ GHz），$f/f_p\in\{0.5,0.9,1.1,2\}$；
- 8a：PIC（$\ge16$ 粒子/格，形状 2 阶；$\Delta x\le\lambda_D/2$ 检查）；
- 8b：冷流体（`fluids.species_names`，负电荷电子流体）；
- 测量：两点相位差 $\Rightarrow k(\omega)$；截止附近反射；趋肤深度。
- 8c：仪器地板（无等离子体）。

**测量方法**：两探针 $E_z$ 的时间序列互谱相位差；拟合 $\omega(k)$；
反射用 C1 的双波分解；电荷/电流守恒用 `ChargeOnEB`/总电荷诊断。

**判据**：$\omega(k)$ 与解析 $\le1.5\%$；截止频率 $\le2\%$；趋肤深度 $\le5\%$；
PIC 与冷流体两条路径互相 $\le3\%$；电荷守恒 $\le10^{-6}$（相对）。
**交付物**：`cases/c08_plasma/`（色散图、截止、两条路径对照）。
**风险**：PIC 噪声 → 用大粒子数 + 时间平均；流体路径的边界条件（出流）会反射 → 用 PML 或长域。

**文献门槛**：Birdsall & Langdon《Plasma Physics via Computer Simulation》；
Stix《Waves in Plasmas》第 1 章；Krall & Trivelpiece。**要复现**：冷等离子体色散曲线
（Stix 图 1-2 一类）的数值表。

---

### C9　磁化等离子体波模：whistler / Alfvén / ion-Bernstein

**目标**：验证含背景磁场的波模，这是 MHD/混合模型的直接物理基础；同时标定 hybrid（Ohm）
与 Darwin 两条求解器路径。

**依赖**：`algo.maxwell_solver = hybrid` 或 `darwin`（现成）。示例 `Examples/Tests/magnetized_plasma_modes/`
提供 PICMI 输入（当前 PYTHON=OFF，需改写为 `inputs` 或打开 PYTHON）。

**物理与理论解**（附录 A.8）
- 平行传播（$\mathbf{k}\parallel\mathbf B_0$）：$n^2=1-\dfrac{\omega_p^2}{\omega(\omega\mp\omega_c)}$，
  $-$ 为 R 模、$+$ 为 L 模；whistler（$\omega\ll\omega_{ce}$）：
  $n^2\approx\dfrac{\omega_p^2}{\omega\omega_{ce}}$ ⇒ $\omega\approx\dfrac{c^2k^2\omega_{ce}}{\omega_p^2}$。
- Alfvén（$\omega\ll\omega_{ci}$）：$\omega=kv_A\cos\theta$，$v_A=\dfrac{B_0}{\sqrt{\mu_0\rho}}$。
- ion-Bernstein：$\omega$ 接近 $\omega_{ci}$ 的高次谐波（Stix）。
- 本算例同时检查：$\mathbf E$、$\mathbf B$ 的极化比（R/L 模的圆极化）。

**数值设置**
- 1D/2D：$B_0$ 沿 $z$ 或 $x$；密度 $10^{18}$–$10^{19}\ \mathrm{m^{-3}}$；$T_e=T_i=100$ eV；
  $\omega/\omega_{ci}\in\{0.1,0.3,1,3\}$；测 $\omega(k)$ 与极化。
- 9a hybrid；9b Darwin；9c 仪器地板（$B_0=0$，退化为 C8）。
- 时间 ≥ 200 $\omega_{ci}^{-1}$。

**测量方法**：时空 FFT（$k$ 由空间 FFT，$\omega$ 由时间 FFT）；峰位拟合；
极化比 $|E_x|/|E_y|$ 与相位差。
**判据**：峰位 $\le2\%$；极化比 $\le5\%$；hybrid 与 Darwin 互差 $\le5\%$；无色散泄漏。
**交付物**：`cases/c09_magnetized/`（色散图、极化、两条求解器对照）。
**风险**：hybrid 是 Darwin/Ohm（无位移电流）→ 高频模（$\omega\gtrsim\omega_{pe}$）不可信，
必须在报告中画出"模型适用边界"。**这条要写进结论，不能当失败隐藏**。

**文献门槛**：Stix《Waves in Plasmas》第 1 章（冷等离子体波）；WarpX `magnetized_plasma_modes` 示例的
`README.rst` 与 `analysis.py`（项目内既有实现，作为第二套对照）；Wang et al. 2015
（DOI 10.1063/1.4906063，5 矩/10 矩 vs PIC 的基准，供 C10 用）。

---

### C10　核壳球（Aden–Kerker）+ 曲面 cut-cell + AMR：综合验收

**目标**：把前 9 个算例的自由度**同时**放进来：曲面、两层材料、窄共振、AMR、长时间守恒。
这是"能不能用于真实介质/等离子体问题"的最终门槛。

**依赖**：M2（介质 + EB 的场更新）与 M4（每层材料 + AMR）。**现在不能跑**。

**物理与理论解**（附录 A.9）
- 核（$n_1$，半径 $a_1$）+ 壳（$n_2$，外半径 $a_2$）+ 背景 $n_3$：Aden–Kerker 双层 Mie 系数
  （递推形式），给出 $Q_{\rm ext}(x)$、$S(\theta)$、近场。
- 两个严格极限用于交叉验证：
  1. 壳去掉（$n_2=n_3$）⇒ 单球 Mie（C4）；
  2. 核去掉（$n_1=n_2$）⇒ 均匀球 Mie（C4）。
- 共振：壳厚变化时 $Q_{\rm ext}$ 的峰位移动；AMR 收敛要求"细层加密不改变峰位超过 0.2%"。

**数值设置**
- 10a（无 AMR）：$x=k a_2\in\{1,2,5\}$（即 $a_2=x\lambda_0/2\pi$），$a_1=0.5a_2$，$n_1=1.5$，$n_2=2.0$。
- 10b（AMR）：在界面附近加密 1–2 层（refinement ratio 2），比较粗细层结果；
  `amr.max_level={0,1,2}`；要求 AMR 不引入伪反射（对比 10a 的远场）。
- 10c（极限校验）：$n_2=n_3$ 与 $n_1=n_2$ 两个退化配置，必须回到 C4 的结果（$\le2\%$）。
- 10d（可选，磁化等离子体背景）：把背景换成 C9 的磁化等离子体；无解析解，
  只用**守恒 + 两个极限**验收（等离子体关掉→10a；壳关掉→C9）。

**测量方法**：同 C4（闭合面坡印廷通量 + 近场 L2）；AMR 情形另做 coarse–fine 通量守恒；
能量守恒用 `FieldEnergy`，电荷守恒用总电荷与 `ChargeOnEB`。

**判据**：$Q_{\rm ext}\le2\%$；近场 $\le3\%$；极限校验 $\le2\%$；
AMR 峰位漂移 $\le0.2\%$；能量漂移 $\le10^{-3}$（100 周期）；电荷守恒 $\le10^{-6}$。
**交付物**：`cases/c10_coated_sphere/`（含 AMR 收敛表、极限校验、守恒曲线）。
**风险**：M2/M4 未完成 → 本算例保持红色，直到 §8.3 的 P1/P2 通过。

**文献门槛**：Aden & Kerker 1951（DOI 10.1063/1.1699834）；Bohren & Huffman
（核壳球级数）；Wang et al. 2015（多流体矩 vs PIC，供 10d 的物理解释）。

---

## 5. 依赖图与执行顺序

```
            ┌──────────────┐
            │ C1 Fresnel   │  材料 + 斜入射几何（倾斜界面）
            └──────┬───────┘
                   │
        ┌──────────┴───────────┐
        ▼                      ▼
  ┌────────────┐         ┌────────────┐
  │ C2 cut-cell│         │ C3 Fabry–P │
  │ (介质, M2) │         └─────┬──────┘
  └─────┬──────┘               │
        │                      ▼
        │              ┌──────────────┐
        │              │ C5 lossy     │
        │              └──────┬───────┘
        │                     ▼
        │              ┌──────────────┐
        │              │ C6 multilayer│ (薄层 → M2)
        │              └──────┬───────┘
        ▼                     ▼
  ┌───────────┐        ┌──────────────┐
  │ C4 Mie    │◄───────┤ C10 coated   │ (M2 + M4)
  │ PEC→介质  │        │ + AMR        │
  └─────┬─────┘        └──────┬───────┘
        │                     ▲
        ▼                     │
  ┌───────────┐        ┌──────────────┐
  │ C7 TIR    │        │ C8/C9 plasma │
  └───────────┘        └──────────────┘
```

- **第一波（不依赖新功能）**：C1（1a/1b/1c/1d/1e）→ **C2-0**（ECT 几何自检）→ C3 → C5 → C7 → C8 → C9。
  C2 本体（介质 cut-cell）在 M2 之后，与 C1 成对运行。
  这一波可以在本周内产出可复查的基线。
- **第二波（依赖 M2）**：C2（与 C1 逐参数配对）、C4-介质、C6-薄层、C10a/10b/10c。
- **第三波（依赖 M4 + 耦合）**：C10d。
- 每个算例的"仪器地板"子算例与主算例**同一天**跑，否则本底不可比。

---

## 6. 交付物规范

### 6.1 目录树

```
benchmarks/
  STATUS.md                 # 矩阵：C1–C10 × {未跑/通过/失败/阻塞}
  CHANGELOG.md              # 阈值/方法/结论的每次变更
  tools/
    report.py               # 从各 case 的 results.json 生成总报告
    verify.py               # 从 raw/ 重算并比对 results.json
    test_theory.py          # 解析解单测（pytest）
  cases/
    c01_fresnel/
      inputs                # 实际运行输入（唯一真源）
      gen.py                # 从参数表生成 inputs（含随机种子固定）
      theory.py             # 解析解（与 analyze.py 分离）
      analyze.py            # 从 raw/ 计算 r,t,θ_t,能量
      expected.json         # 预注册阈值
      run.sh                # 运行命令（含 env -u LD_LIBRARY_PATH）
      meta.json             # 运行溯源
      sha256sums.txt        # raw/ 与二进制哈希
      raw/                  # plotfile / reduced diags（只读）
      results.json          # analyze.py 输出（唯一数据源）
      figs/                 # 由 results.json 生成
      LITERATURE.md         # G0 文献卡
      README.md             # 结论摘要（不写数值，数值在 results.json）
```

### 6.2 `meta.json` 必填字段

```json
{
  "case": "c01_fresnel",
  "timestamp_utc": "…",
  "hostname": "…",
  "binary": {"path": ".../warpx.3d...", "sha256": "…"},
  "inputs_sha256": "…",
  "cmake_cache": {"WarpX_DIMS": "3", "WarpX_EB": "ON", "WarpX_MPI": "ON",
                  "WarpX_PYTHON": "OFF", "WarpX_OPENPMD": "OFF"},
  "compiler": "gcc --version 输出首行",
  "mpi": "mpirun --version 输出首行",
  "command": "mpiexec -n 8 ./warpx.3d inputs",
  "walltime_s": 0,
  "source_sha256_manifest": "sha256sums_source.txt"
}
```

> 注意：`WarpX_fluid` 当前 git **无提交**，不能用 commit 溯源；必须用 `Source/` 全量 sha256 清单。

### 6.3 `results.json` 必填字段（示例）

```json
{
  "case": "c01_fresnel",
  "subcases": [
    {"id": "1a_eps4_normal", "eps_r": 4.0, "sigma": 0.0,
     "r": {"re": 0.0, "im": 0.0}, "t": {"re": 0.0, "im": 0.0},
     "R": 0.0, "T": 0.0, "r_theory": {"re": 0.0, "im": 0.0},
     "err_r_rel": 0.0, "err_phase_deg": 0.0,
     "window": [0.0, 0.0], "probes": []}
  ],
  "floor": {"case": "1e_vacuum", "err_r_rel": 0.0},
  "convergence": [{"dx": 0.0, "err": 0.0}],
  "verdict": "pass|fail|blocked",
  "notes": "…"
}
```

所有浮点数保留**原始精度**（不四舍五入），图与报告各自格式化。

### 6.4 报告生成规则

- `tools/report.py` 只读 `results.json` + `expected.json`；
- 生成 `REPORT.md`：每个算例一节（判据、实测、差值、结论）；
- **禁止**在 `REPORT.md` 里出现任何手工数字；发现手工数字即视为数据污染；
- 图由 `results.json` 直接生成，图注自动写入数据来源（case id + timestamp）。

---

## 7. 汇总判据表

| 算例 | 主要物理量 | A 级 | B 级 | 本底要求 |
|---|---|---|---|---|
| C1 | $\lvert r\rvert,\lvert t\rvert,\theta_t,R+T$ | 2%, 0.8°, 1% | 1%, 0.5°, 0.5% | $<10^{-4}$ |
| C2 | 对 Fresnel 的 r/t、与 C1 之差、偏移敏感性、收敛 | 2%, 1%, 1%, ≥1.5 | 1%, 0.5%, 0.5%, ≥1.8 | $<10^{-4}$ |
| C2-0 | $\lvert r+1\rvert$、$\epsilon_{\rm spur}$、能量 | 1%, $10^{-3}$, $10^{-3}$ | 0.5%, $10^{-4}$, $10^{-4}$ | 同 |
| C3 | $T(f)$、峰位、FSR、finesse | 2%, 0.5%, 1%, 5% | 1%, 0.3%, 0.5%, 3% | $<10^{-4}$ |
| C4 | $Q_{\rm ext}$、近场、$S(\theta)$ | 3%, 5%, 8% | 2%, 3%, 5% | $Q<10^{-3}$ |
| C5 | $\delta_s$、$R,T,A$ | 3%, 2%, 1% | 2%, 1%, 0.5% | $<10^{-4}$ |
| C6 | $R(f)$、$R_{\min}$、阻带中心 | 2%, 0.5%, 1% | 1%, 0.3%, 0.5% | $<10^{-4}$ |
| C7 | $R$、$\kappa$、PML 伪反射 | 0.1%, 3%, $10^{-3}$ | 0.05%, 2%, $10^{-4}$ | 同 |
| C8 | $\omega(k)$、截止、$\delta_s$ | 1.5%, 2%, 5% | 1%, 1%, 3% | 无等离子体 |
| C9 | 峰位、极化、两求解器差 | 2%, 5%, 5% | 1%, 3%, 3% | $B_0=0$ |
| C10 | $Q_{\rm ext}$、近场、极限、守恒 | 2%, 3%, 2%, $10^{-3}$ | 1%, 2%, 1%, $10^{-4}$ | — |

> 阈值是**预注册初值**；若 G0/G1 发现文献本身给出的不确定度更大（例如 Mie 在 $x<1$ 的收敛慢），
> 允许通过 `CHANGELOG.md` 调整，但必须写明理由与文献依据，并在报告中同时给出新旧阈值下的结论。

---

## 8. 风险、回退与能力缺口

### 8.1 已识别的数值风险

| 风险 | 影响算例 | 对策 |
|---|---|---|
| 阶梯 vs cut-cell 混用导致误差来源不清 | C1–C7 | 每个算例都同时跑"阶梯对照"，明确当前路径 |
| PML 对倏逝波吸收不完全 | C7 | PML 厚度/距离扫描；必要时用"长真空 + 大域"替代 |
| 未填充域被误判为不稳定 | 全部 | 只取波前通过后的时间窗；用 $\max$ 而非全局能量 |
| 2D 圆柱 $Q_{\rm ext}$ 归一化约定 | C4 | G1 用 Bohren & Huffman 表逐点核对后固化 |
| hybrid/Darwin 无位移电流 | C9 | 报告里画"模型适用边界"，不把高频失效当 bug |
| 1D/2D 材料路径未验证 | C1 | 1a 先跑 1D；失败则退 3D 薄板 |
| 曲面 ECT 最小面分数过小 | C2/C4 | 记录 `area_mod` 统计；必要时局部加密（C10） |

### 8.2 数据完整性风险

| 风险 | 对策 |
|---|---|
| 手工改数/复制粘贴 | `verify.py` 重算 + 哈希清单 + 报告自动生成 |
| 事后调阈值 | `expected.json` 预注册 + `CHANGELOG` |
| 删除失败算例 | `STATUS.md` 强制登记；`raw/` 只读 |
| 环境漂移（库版本） | `meta.json` 记录编译器/MPI；关键库版本写入 |
| 无 commit 可溯源 | `Source/` 全量 sha256 清单 |

### 8.3 能力缺口与前置改造（对应旧调研的 M2/M4）

| 编号 | 缺口 | 影响 | 验收（进入本计划的条件） |
|---|---|---|---|
| P1（M2） | 介质 + EB 的场更新（EB 内继续解 Maxwell、界面平均） | **C2**、C4-介质、C6-薄层、C10 | C2 与 C1 同参数之差 ≤1%；在 C4 的 PEC 基线上加介质，$Q_{\rm ext}$ 达 C4 判据；C3 亚网格 $T(f)$ 达 2% |
| P2（M4） | 每层材料系数 + AMR 材料界面 | C10b | 10a 与 10b 峰位差 $\le0.2\%$；无 coarse–fine 伪反射 |
| P3 | 斜入射源（TF/SF 或相位斜坡 laser） | 可选（本计划用倾斜界面替代） | 与倾斜界面结果互差 $\le1\%$ |
| P4 | 近远场变换（NTFF） | 可选（本计划用近场闭合面） | 与近场法互差 $\le2\%$ |

> **重要**：P1 是"介质 + cut-cell"的唯一阻塞项。按旧调研文档 §3.4 的结论，正确路线是
> "材料平均的 Yee（Meep 式 subpixel）+ EB 几何感知"，而**不是**把 ECT 的 PEC 公式改成介质版；
> 后者会同时破坏两种语义。本计划的 C2-0（PEC/ECT 极限）与 C2/C4-介质（材料 Yee）分开，就是为了
> 让两条语义各自可验证。

---

## 9. 文献（均已在 G0 中核对 DOI）

**电磁/数值**
1. S. Dey, R. Mittra, *A locally conformal finite-difference time-domain (FDTD) algorithm for
   modeling three-dimensional perfectly conducting objects*, IEEE Microwave and Guided Wave
   Letters 7(9), 273–275 (1997). DOI: [10.1109/75.622536](https://doi.org/10.1109/75.622536)
2. T. Xiao, Q. H. Liu, *Enlarged cells for the conformal FDTD method to avoid the time step
   reduction*, IEEE Microwave and Wireless Components Letters 14(12), 551–553 (2004).
   DOI: [10.1109/LMWC.2004.837384](https://doi.org/10.1109/LMWC.2004.837384)
3. A. Farjadpour *et al.*, *Improving accuracy by subpixel smoothing in the finite-difference
   time domain*, Optics Letters 31(20), 2972–2974 (2006).
   DOI: [10.1364/OL.31.002972](https://doi.org/10.1364/OL.31.002972)
4. A. Deinega, I. Valuev, *Subpixel smoothing for conductive and dispersive media in the
   finite-difference time-domain method*, Optics Letters 32(23), 3429–3431 (2007).
   DOI: [10.1364/OL.32.003429](https://doi.org/10.1364/OL.32.003429)
5. K. R. Umashankar, A. Taflove, *A novel method to analyze electromagnetic scattering of
   complex objects*, IEEE Trans. Electromagnetic Compatibility 24(4), 397–405 (1982).
   DOI: [10.1109/TEMC.1982.304054](https://doi.org/10.1109/TEMC.1982.304054)

**散射/光学**
6. G. Mie, *Beiträge zur Optik trüber Medien, speziell kolloidaler Metallösungen*,
   Annalen der Physik 330(3), 377–445 (1908).
   DOI: [10.1002/andp.19083300302](https://doi.org/10.1002/andp.19083300302)
7. C. F. Bohren, D. R. Huffman, *Absorption and Scattering of Light by Small Particles*,
   Wiley (1983/1998). DOI: [10.1002/9783527618156](https://doi.org/10.1002/9783527618156)
8. W. J. Wiscombe, *Improved Mie scattering algorithms*, Applied Optics 19(9), 1505–1509 (1980).
   DOI: [10.1364/AO.19.001505](https://doi.org/10.1364/AO.19.001505)
9. A. L. Aden, M. Kerker, *Scattering of electromagnetic waves from two concentric spheres*,
   Journal of Applied Physics 22(10), 1242–1246 (1951).
   DOI: [10.1063/1.1699834](https://doi.org/10.1063/1.1699834)
10. H. A. MacLeod, *Thin-Film Optical Filters*, 4th ed., CRC Press (2010).
    DOI: [10.1201/9781420073034](https://doi.org/10.1201/9781420073034)

**等离子体/矩方程**
11. L. Wang, A. H. Hakim, A. Bhattacharjee, K. Germaschewski, *Comparison of multi-fluid moment
    models with particle-in-cell simulations of collisionless magnetic reconnection*,
    Physics of Plasmas 22, 012108 (2015).
    DOI: [10.1063/1.4906063](https://doi.org/10.1063/1.4906063)
12. S. I. Braginskii, *Transport processes in a plasma*, Reviews of Plasma Physics 1 (1965)
    （碰撞/黏性闭包；如用到）。
13. N. A. Krall, A. W. Trivelpiece, *Principles of Plasma Physics*（色散关系标准参考）。
14. T. H. Stix, *Waves in Plasmas*, AIP (1992)（磁化等离子体色散）。
15. C. K. Birdsall, A. B. Langdon, *Plasma Physics via Computer Simulation*（PIC 噪声与色散）。

> 说明：12–15 为书籍/综述，DOI 不适用；G0 中直接标注版次与章节。
> 所有 DOI 均通过 Crossref 检索核对（2026-10-31），链接可直接解析。

---

## 附录 A　理论解（供 theory.py 实现）

### A.1 复 Fresnel（$\mu_r=1$）

设 $\tilde\varepsilon_j=\varepsilon_{rj}\varepsilon_0-i\sigma_j/\omega$，
$\tilde n_j=\sqrt{\tilde\varepsilon_j/\varepsilon_0}$，
$k_0=\omega/c$，$k_x=k_0\tilde n_1\sin\theta_i$，
$k_{zj}=k_0\sqrt{\tilde n_j^2-\tilde n_1^2\sin^2\theta_i}$（分支：$\mathrm{Re}\ge0$，$\mathrm{Im}\ge0$）。

$$
r_s=\frac{k_{z1}-k_{z2}}{k_{z1}+k_{z2}},\quad
t_s=\frac{2k_{z1}}{k_{z1}+k_{z2}},
$$
$$
r_p=\frac{\tilde\varepsilon_2k_{z1}-\tilde\varepsilon_1k_{z2}}{\tilde\varepsilon_2k_{z1}+\tilde\varepsilon_1k_{z2}},\quad
t_p=\frac{2\tilde n_1\cos\theta_i}{\tilde n_2\cos\theta_i+\tilde n_1\cos\theta_t}.
$$

功率：$R_{s,p}=|r_{s,p}|^2$，$T_{s,p}=\dfrac{\mathrm{Re}\,k_{z2}}{\mathrm{Re}\,k_{z1}}|t_{s,p}|^2$。

### A.2 单层板（Airy）

$$
\delta=k_{z2}d=k_0\tilde n_2d\cos\theta_t,\quad
r=\frac{r_{12}+r_{23}e^{2i\delta}}{1+r_{12}r_{23}e^{2i\delta}},\quad
t=\frac{t_{12}t_{23}e^{i\delta}}{1+r_{12}r_{23}e^{2i\delta}}.
$$

无耗、$n_1=n_3$：$r_{23}=-r_{12}$，$T=\dfrac{(1-R)^2}{(1-R)^2+4R\sin^2\delta}$。

### A.3 多层特征矩阵

$$
M_j=\begin{pmatrix}\cos\delta_j & \dfrac{i}{\eta_j}\sin\delta_j\\[4pt]
i\eta_j\sin\delta_j & \cos\delta_j\end{pmatrix},\quad
\delta_j=k_0n_jd_j\cos\theta_j,\quad
\eta_j=\begin{cases}n_j\cos\theta_j & s\\ n_j/\cos\theta_j & p\end{cases}
$$
$$
\begin{pmatrix}B\\C\end{pmatrix}=\left(\prod_{j=N}^{1}M_j\right)\begin{pmatrix}1\\\eta_s\end{pmatrix},\quad
r=\frac{\eta_0B-C}{\eta_0B+C},\quad t=\frac{2\eta_0}{\eta_0B+C},\quad T=\frac{\eta_s}{\eta_0}|t|^2.
$$

### A.4 球 Mie

$$
\psi_n(x)=xj_n(x),\quad \xi_n(x)=xh_n^{(1)}(x),\quad m=\tilde n_s/\tilde n_m,\quad x=k a.
$$
$$
a_n=\frac{m\psi_n(mx)\psi_n'(x)-\psi_n(x)\psi_n'(mx)}
          {m\psi_n(mx)\xi_n'(x)-\xi_n(x)\psi_n'(mx)},\quad
b_n=\frac{\psi_n(mx)\psi_n'(x)-m\psi_n(x)\psi_n'(mx)}
          {\psi_n(mx)\xi_n'(x)-m\xi_n(x)\psi_n'(mx)}.
$$
$$
Q_{\rm ext}=\frac{2}{x^2}\sum_{n=1}^{\infty}(2n+1)\mathrm{Re}(a_n+b_n),\quad
Q_{\rm sca}=\frac{2}{x^2}\sum_{n=1}^{\infty}(2n+1)(|a_n|^2+|b_n|^2),\quad
Q_{\rm abs}=Q_{\rm ext}-Q_{\rm sca}.
$$
PEC 极限：$a_n=\psi_n'(x)/\xi_n'(x)$，$b_n=\psi_n(x)/\xi_n(x)$。

### A.5 圆柱（正入射）

$x=ka$，$m=\tilde n_s/\tilde n_m$，$J_n$ 为 Bessel、$H_n^{(1)}$ 为 Hankel：
$$
a_n=\frac{mJ_n(mx)J_n'(x)-J_n(x)J_n'(mx)}
          {mJ_n(mx)H_n^{(1)\prime}(x)-H_n^{(1)}(x)J_n'(mx)}\quad(\text{TM}, E\parallel\text{轴}),
$$
$$
b_n=\frac{J_n(mx)J_n'(x)-mJ_n(x)J_n'(mx)}
          {J_n(mx)H_n^{(1)\prime}(x)-mH_n^{(1)}(x)J_n'(mx)}\quad(\text{TE}, H\parallel\text{轴}).
$$
2D 截面（单位长度）：$C_{\rm ext}^{\rm TM}=\dfrac{4}{k}\sum_{n=-\infty}^{\infty}\mathrm{Re}\,a_n$，
$Q_{\rm ext}=C_{\rm ext}/(2a)$；TE 同理用 $b_n$。
> **归一化常数在 G1 中与 Bohren & Huffman 表逐点核对后固化**；不得凭记忆写死。

### A.6 全内反射 / 有耗

$$
\sin\theta_c=\frac{n_2}{n_1}\ (n_1>n_2),\qquad
\kappa=k_0\sqrt{n_1^2\sin^2\theta_i-n_2^2},\qquad E\propto e^{-\kappa z}.
$$
有耗：$\tilde n=\sqrt{\varepsilon_r(1-i\tan\delta)}$，$\kappa=\mathrm{Im}(k_0\tilde n)$，
$\delta_s=1/\kappa$；良导体 $\delta_s=\sqrt{2/(\omega\mu\sigma)}$。

### A.7 无磁化等离子体

$$
\omega_p=\sqrt{\frac{n_ee^2}{\varepsilon_0m_e}},\qquad
\omega^2=\omega_p^2+c^2k^2,\qquad
\kappa=\frac{\sqrt{\omega_p^2-\omega^2}}{c}\ (\omega<\omega_p),\qquad
r=\frac{k_0-k}{k_0+k}.
$$

### A.8 磁化等离子体（Stix，平行传播）

$$
n^2=1-\frac{\omega_p^2}{\omega(\omega\mp\omega_c)}\quad(R/L),\qquad
\omega_{ce}=\frac{eB_0}{m_e},\quad \omega_{ci}=\frac{ZeB_0}{m_i}.
$$
whistler（$\omega\ll\omega_{ce}$）：$n^2\approx\dfrac{\omega_p^2}{\omega\omega_{ce}}$
⇒ $\omega\approx\dfrac{c^2k^2\omega_{ce}}{\omega_p^2}$。
Alfvén：$v_A=\dfrac{B_0}{\sqrt{\mu_0\rho}}$，$\omega=kv_A\cos\theta$。

### A.9 核壳球（Aden–Kerker）

双层球系数按 Aden–Kerker 递推：核层用 $(m_1x_1)$ 的 Riccati–Bessel 解，壳层在其上做
第二层匹配，得到外层 $a_n,b_n$；$Q_{\rm ext}$ 仍用 A.4 的求和式。
**实现要求**：`theory.py` 先复现 Aden–Kerker 文中 $m_1=1.5,m_2=2.0$ 的 $Q_{\rm ext}(x)$ 表，
再用于 C10；极限 $n_2=n_3$ 与 $n_1=n_2$ 必须自动退化为 A.4。

---

## 附录 B　输入模板

### B.1 C1a：1D 正入射（材料）

```
geometry.dims = 1
geometry.prob_lo = -4.0
geometry.prob_hi =  4.0
amr.n_cell = 1600
amr.max_level = 0
boundary.field_lo = pml
boundary.field_hi = pml
warpx.cfl = 0.5
algo.maxwell_solver = yee
algo.em_solver_medium = macroscopic
macroscopic.epsilon = 4.0*epsilon0
macroscopic.mu = mu0
macroscopic.sigma = 0
lasers.names = src
src.position = 0.0
src.direction = 1.0
src.polarization = 0.0 1.0 0.0
src.wavelength = 1.0
src.e_max = 1.0e6
src.profile = parse_field_function
src.field_function(X,Y,t) = "exp(-((t-t0)/tau)**2)*cos(omega0*t)"
diagnostics.diags_names = d1
d1.intervals = 1
d1.diag_type = Full
d1.fields_to_plot = Ey Bz
```
> 实际 `inputs` 由 `gen.py` 生成，界面位置、$\varepsilon_r$、$\sigma$、探针位置全部参数化。

### B.2 C2：2D 倾斜界面（介质 cut-cell；M2 前只能用 C2-0）

**B.2a　C2-0（PEC 极限，现在可跑）**

```
geometry.dims = 2
geometry.prob_lo = -10.0 -10.0
geometry.prob_hi =  30.0  30.0
amr.n_cell = 1600 1600
amr.max_level = 0
boundary.field_lo = pml pml
boundary.field_hi = pml pml
warpx.cfl = 0.4
algo.maxwell_solver = ect
warpx.eb_implicit_function = "cos(alpha)*(x-x0)+sin(alpha)*(z-z0)"   # 倾斜平面
```

**B.2b　C2（介质 cut-cell，M2 之后）**：与 C1 相同的源/域/监测面，界面改为 EB，材料在
EB 内部由 `epsilon_inside_eb`（或亚网格等级 set 采样）给出；**在 M2 落地前此输入不可运行**。
模板字段：
```
algo.em_solver_medium = macroscopic
macroscopic.epsilon_function(x,y,z) = "..."     # 流体侧的 ε
# 介质侧：由 EB 内部材料提供（M2 的接口名在实现时确定，例如 epsilon_inside_eb）
warpx.eb_implicit_function = "cos(alpha)*(x-x0)+sin(alpha)*(z-z0)"
```

> 符号约定：`eb_implicit_function` 的符号决定哪一侧是"物体内部"，必须与材料侧一致；
> 在 `gen.py` 中由 $\alpha$ 与目标侧唯一确定，并写入 `expected.json`。

### B.3 C3：1D 平板

```
algo.em_solver_medium = macroscopic
macroscopic.epsilon_function(x,y,z) = "epsilon0*(1 + (er-1)*(z>0)*(z<d))"
```
> 亚网格子例同一模板，仅把 `d` 改为 $\lambda/20$。

---

## 附录 C　脚本接口约定

### C.1 `theory.py`

```python
def fresnel(n1, n2, theta_i_deg, pol="s") -> dict:
    """返回 {'r': complex, 't': complex, 'R': float, 'T': float,
             'theta_t_deg': float}；n1,n2 为复折射率。"""

def airy(n1, n2, n3, d, lam, theta_i_deg, pol="s") -> dict: ...
def transfer_matrix(layers, lam, theta_i_deg, pol="s") -> dict: ...
def mie_sphere(a, m, x, n_max=...) -> dict: ...
def mie_cylinder(a, m, x, pol="TM", n_max=...) -> dict: ...
def mie_coated(a_core, a_shell, m_core, m_shell, m_bg, x) -> dict: ...
def tir_decay(n1, n2, theta_i_deg, k0) -> float: ...
def plasma_dispersion(ne, f) -> dict: ...
def stix_parallel(ne, B0, f, species="e") -> dict: ...
```

### C.2 `analyze.py`

```python
def load_fields(raw_dir): ...
def two_wave_fit(E, z, t, omega, k): ...   # 返回 A_plus, A_minus
def extract_r_t(ds, cfg): ...
def energy_balance(ds, cfg): ...
def write_results(results: dict, path): ...  # 只写 results.json
```

### C.3 `verify.py`

```text
1. 读 meta.json，核对二进制/输入哈希；
2. 重跑 analyze.py 到临时目录；
3. 逐字段比较 results.json（相对误差 <1e-12）；
4. 不一致 → exit 1。
```

---

## 附录 D　算例预算（估算，非承诺）

| 算例 | 典型运行数 | 单次规模 | 备注 |
|---|---|---|---|
| C1 | 5 子例 × 3 分辨率 ≈ 40 | 1D/2D，分钟级 | 1D 最快 |
| C2 | 8 角度 × 3 ε × 2 偏振 ≈ 48，加 σ/偏移/收敛 ≈ 80 | 2D 1600²，~10 分钟/次 | 每例与 C1 配对 |
| C2-0 | 8 角度 × 3 分辨率 ≈ 24 | 2D，ECT 单层 | 现在可跑 |
| C3 | 3 厚度 × 3 分辨率 + 斜入射 ≈ 20 | 1D/2D | |
| C4 | 5 $x$ × 2 偏振 × 3 分辨率 × 2（PEC/介质）≈ 60 | 2D/3D，小时级 | 3D 球最贵 |
| C5 | 4 σ × 3 分辨率 ≈ 15 | 1D | 高 σ 需很细网格 |
| C6 | 2 设计 × 3 分辨率 × 2 角度 ≈ 15 | 1D | |
| C7 | 4 角度 × 3 分辨率 ≈ 15 | 2D | |
| C8 | 4 频率 × 2 路径 × 2 分辨率 ≈ 16 | 1D/2D | |
| C9 | 4 频率 × 2 求解器 × 2 分辨率 ≈ 16 | 1D/2D | |
| C10 | 3 $x$ × 3 AMR 级别 × 3 分辨率 ≈ 27 | 3D，最贵 | M2/M4 后 |

总计约 250–300 次运行；在单机 8 核 OMP 上按"分钟–小时"量级估计，需要分批执行与断点续跑
（`run.sh` 支持 `--resume`）。

---

## 附录 E　立即可以做的三件事（不依赖 M2/M4）

1. **建 1D 与 2D 构建**：`cmake -S . -B build_1d -DWarpX_DIMS=1 -DWarpX_EB=ON -DWarpX_MPI=ON`；
   `build_2d` 同理（`WarpX_DIMS=2`）。先跑 C1a 的 1D 正入射，确认材料路径在 1D 可用。
2. **把 `magnetized_plasma_modes` 的 PICMI 输入翻译成 `inputs`**（当前 `WarpX_PYTHON=OFF`），
   作为 C9 的起点。
3. **写 `theory.py` + `test_theory.py`**：先实现 A.1（Fresnel）与 A.4（Mie），
   用 §9 的文献表做 G1 校验；这一步与仿真无关，可以并行做。
