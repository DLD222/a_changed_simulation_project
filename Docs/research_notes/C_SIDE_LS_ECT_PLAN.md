# C 侧 G1/G2 实施方案：level-set 精确 `Lx/Lz` 与 TE 面完整 ECT

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 日期：2026-10-31 · 前置：`C_SIDE_FINDINGS_2026-10-31.md`（F-1/F-2/F-3）、`C_SIDE_FIX_PLAN.md`、
> `C_SIDE_ORDER_ANALYSIS_REVIEW.md`、`C_SIDE_DESIGN2.md` §1.5（方案 C）
> 范围：本报告只讲两件事——**① `Lx/Lz` 的 level-set 精确积分怎么实现（G1）；② TE 面（`Bx`/`Bz`）
> 怎么启用完整 ECT 借用（G2）**，以及配套的稳定化（G3）与验收顺序。不改代码，行号指当前 `master`。
>
> 结论速览：
> - G1 的核心不是"换公式"，而是**把 `Lx/Lz` 从"2D 面积分数估 1D 线长"改成"沿对偶线与 EB 求交的 1D 积分"**，
>   维度匹配后几何误差从 O(h) 变成 O(h²)（用解析 φ 则近机器零）；
> - G2 不是只修 `EvolveB` 索引一处，而是 **4 处配套 + 1 处索引**；邻居编码约定必须在写代码前定死并用
>   单 cell 真值表验证；
> - G3（完整稳定化）的正确形态是"分母下限（三个量）+ 阈值冻结（DM 式）+ 掩膜统一 `A_y<=A_stab`"，
>   而不是把 `A_stab` 抬大。

---

## 0. 现状（改动前的锚点）

| 项 | 现状 | 位置 |
|---|---|---|
| `Lx/Lz` 构造 | 方案 B：`½h(v+v)`，一阶（维度错配） | `EmbeddedBoundaryInit.cpp:798-812` |
| `A_y` | `¼h²Σv`，二阶，不动 | 同上 |
| TE 面（`Bx/Bz`）稳定化 | BCK-only，`flag_ext_face=0`，无借用 | `EmbeddedBoundaryInit.cpp:465-482` |
| 扩展循环（XZ） | `idim=1..`（By 才扩展） | `WarpXFaceExtensions.cpp:671`、`:831` |
| 借用邻居映射（`EvolveB`） | idim=0/2 用 3D-y 索引 | `EvolveB.cpp:305/314/337/346` |
| `GetNeigh/SetNeigh/GetNeighPtr` | dim=0 为 3D 写法 | `WarpXFaceExtensions.cpp:277-408` |
| BCK 兜底（XZ） | 0/2 被 `#if` 禁用 | `WarpXFaceExtensions.cpp:612/628` |
| E 侧稳定化 | 只有 `A_y` 分母下限 `A_stab`（默认 0.1·dx·dz） | `EmbeddedBoundaryInit.cpp:825-831` |

---

## 1. G1：level-set 精确 `Lx/Lz`（方案 C）

### 1.1 原理

```
Lx(i, j+½) = ∫_{x_i−h/2}^{x_i+h/2} [φ(x, z_{j+½}) ≤ 0] dx        // 沿 x 的开放长度（Bx 处）
Lz(i+½, j) = ∫_{z_j−h/2}^{z_j+h/2} [φ(x_{i+½}, z) ≤ 0] dz        // 沿 z 的开放长度（Bz 处）
```

- **符号**：WarpX 约定 `eb_implicit_function` **正值为固体** → 开放 = `φ ≤ 0`（`C_SIDE_RESEARCH2.md` §2.3）。
- **为什么能到二阶**：这是 1D↔1D 的求交积分，不再用 2D 面积分数；`A_y` 保持 `¼h²Σv`（已是二阶）。

### 1.2 数据源与"同源"选择（E7 的关键决策）

| 方案 | 数据 | 精度 vs 真几何 | 与 B 侧的一致性 |
|---|---|---|---|
| **A（推荐）** | `eb_fact.getLevelSet()`，在每条线路上做**线性插值**（等价于 AMReX EB2 每 cell 的线性化 level set） | O(h²)（EB2 每 cell 平面近似真圆） | **同源**（与 `face_areas`/`edge_lengths` 来自同一剖分） |
| B | 直接调用 `eb_implicit_function` 的 parser 逐点求值 | 对解析曲面近精确 | **不同源**：与 AMReX 每 cell 平面差 O(h²) |

- `getLevelSet()`：AMReX 工厂暴露（`WarpX.H:77` → `m_ebdc->getLevelSet()`），是**全网格 nodal MultiFab**
  （`AMReX_EBDataCollection.cpp:31-33`：`convert(a_ba, IntVect::TheUnitVector())`，nGrow=`m_ngrow[0]`），
  规则区也有值 → 可以直接支撑"规则域返回精确 h"。
- **建议先做方案 A**：同源把 E7 风险消掉；若后续发现 O(h²) 平面误差仍是瓶颈，再切 B（并把 B 侧也换成
  parser 同源积分）。
- 注意 `Lx` 在 `z=z_{j+½}`（z 向格心），而 level set 在节点上 → 对每个采样点 x 取
  `φ_c = ½(φ(x, z_j) + φ(x, z_{j+1}))`（z 向线性）；`Lz` 对称地取 x 向线性
  `φ_c = ½(φ(x_i, z) + φ(x_{i+1}, z))`。

### 1.3 算法（只跑一次，初始化时，CPU）

```text
ComputeDualGeometry(scheme == 4) 分支：
  for 每条对偶线：
    1) 若整段 φ_c ≤ 0（规则域）     → L = h（精确返回，保证 L2a 逐位）
    2) 若整段 φ_c > 0（完全覆盖）    → L = 0
    3) 否则：N = 8..16 个采样点，逐点求 φ_c，
       找符号变化区间，线性插值求交点（更稳可用 2–3 次二分），
       L = Σ 开放子段长度
  Lx/Lz 算完后照旧施加 A_stab 下限（0 保持 0）。
```

- 复杂度：每层 O(n_cell × N)，640² 秒级；**不进每步内核**。
- `scheme` 参数已存在（`ComputeDualGeometry(..., int scheme, Real A_stab)`），当前 3/4 在
  `WarpX.cpp:586-596` 被断言挡住；实现 4 后把断言上限放开到 4。
- 实现位置：`EmbeddedBoundaryInit.cpp` 内新增 helper（如 `ComputeDualLengthFromLevelSet(...)`），
  在 `scheme == 4` 时替代方案 B 的 `Lx/Lz` 两段；`A_y` 与 constant/covered 分支不动。

### 1.4 G1 验收

| 判据 | 标准 |
|---|---|
| 纯几何阶数（解析圆，对偶线位置） | `Lx/Lz` RMS 阶 ≥ 1.8（方案 A 预期 2；方案 B 预期近机器零） |
| L2a（规则域） | `max|d|=0` 的等价物：rel ≤ 1e-14，零分量经噪声地板 |
| L3 45° TE | `|r+1|` 与斜率不劣于当前（0.0459 / −1.51%） |
| 几何自检打印 | 规则域 `A_y=dx·dz、Lx=dx、Lz=dz` 精确 |

---

## 2. G2：TE 面启用完整 ECT 借用

### 2.1 要改的 5 处（4 处配套 + 1 处索引）

| # | 文件:行 | 现状 | 改法 |
|---|---|---|---|
| 1 | `WarpXFaceExtensions.cpp:671`、`:831` | one-way / eight-ways 的 XZ 循环 `idim=1..` | 扩到 `idim=0..2`；同时把内层"面内邻居"枚举按 XZ 收敛到：idim=0 只沿 z（偏移 `(0,±1)`）、idim=2 只沿 x（偏移 `(±1,0)`） |
| 2 | `WarpXFaceExtensions.cpp:277-408` | `GetNeigh/SetNeigh/GetNeighPtr` 的 dim=0 是 3D-y 写法 `arr(i, j+i_n, k+j_n)` | XZ 分支：dim=0 → `arr(i, j+j_n, k)`；dim=2 → `arr(i+i_n, j, k)`；dim=1 已有 XZ 分支 |
| 3 | `EmbeddedBoundaryInit.cpp:465-482` | XZ `idim=0/2` 走 BCK-only（`flag_ext_face=0`） | 删除该 XZ 特判，让 0/2 走与 3D 相同的"扩展/available"路径（`ComputeSStab` 的 XZ 分支已在 `WarpXFaceExtensions.cpp:96-133`） |
| 4 | `EvolveB.cpp:305/314/337/346` | idim=0 用 `jp=j+vec(0); kp=k+vec(1)`；else（idim=2）用 `ip=i+vec(0); jp=j+vec(1)` | XZ 分支：idim=0 → `ip=i; jp=j+vec(1); kp=k`；idim=2 → `ip=i+vec(0); jp=j; kp=k`（与 #2 的编码一致） |
| 5 | `WarpXFaceExtensions.cpp:612/628`（可选） | BCK 兜底对 XZ 的 0/2 被 `#if !defined(XZ)` 禁用 | 若要"借用失败 → BCK 兜底"，放开 XZ 禁用 |

### 2.2 邻居编码约定（写代码前必须定死）

- `FaceInfoBox::uint8_to_inds` 给的是**面平面内的 2D 偏移**。
- 对 `idim=0`（Bx 的 x 法向面，per-unit-y 后是 z 段）：邻居只有 z → **存 `(i_n=0, j_n=±1)`**，
  映射 `(i, j+j_n, k)`；对 `idim=2`（Bz）：邻居只有 x → **存 `(i_n=±1, j_n=0)`**，映射 `(i+i_n, j, k)`。
- 该约定必须用 **S1a' 单 cell 真值表**验证：构造一个"S 极小、邻居 S 大"的面，确认 `Venl` 记账与
  `borrowing_size/inds/area` 逐一正确（借用方向错了在 release 下是静默错值）。

### 2.3 风险与验证

- **借用 ≠ 更准**：ECT 借用是"用面积换稳定"，启用后 TE 的精度可能变好也可能持平/变差，必须实测；
- B 侧数值全变：L2b/L3/米氏全部基线重记；先跑一个"小斜角 PEC"对比 `BCK-only vs 借用` 的反射系数与能量；
- 顺序建议：**G1 完成后、G2 单独提交**，不与几何改动混在一起。

---

## 3. G3：完整稳定化（F-3 教训的直接产物）

1. **分母下限覆盖三个量**：`A_y`（已有）、`Lx`、`Lz`（Ex/Ez 的分母，现在裸奔）；阈值保守（≈0.05–0.1·h² / h）；
2. **阈值冻结（Dey–Mittra 式）**：`A_y < A_stab` 的节点直接 `E=0`（PEC 表面切向 E=0，物理正当），
   **不要**把 `A_stab` 抬到 10·dx·dz（F-3：环带 +27%，退化阶梯）；
3. **掩膜统一**：`Ey` 判据从 `ly==0`（太严）改成 `A_y <= A_stab`（比 F3 的 `A_y<=0` 多一个安全边），
   并重开 F3（传回 `dual_A_y`）；
4. E 侧借用（`Venl` 等价物）**暂缓**，等 G1–G3 验证后再评估。

G3 验收（F-2 的教训，三行缺一不可）：

```
max|E|（全网格）                 -- 发现热斑
max|E|（PEC 内 r < 0.9R）        -- 必须为 0 / 机器零
L2（环带 1.5R–2.5R）             -- 原有精度指标
```

---

## 4. 执行顺序与验收总表

| 步 | 内容 | 验收 | 备注 |
|---|---|---|---|
| **G1** | level-set `Lx/Lz`（方案 A 同源版） | 几何阶 ≥1.8；L2a 逐位；L3 不劣化 | 自包含，先做 |
| **G2** | TE 面 ECT 借用（§2 的 5 处） | 单 cell 真值表；L2b/L3 复测；BCK-only vs 借用对比 | 单独提交，基线重记 |
| **G3** | 完整稳定化（分母下限×3 + `A_y<=A_stab` 冻结 + F3 重开） | 三行验收；无热斑；环带不退化 | 依赖 G1 |
| **G4** | 米氏 4 档 × {纯 Yee, ECT+Yee-E, ECT+共形E} + 真空地板 + RCS 数值入射 | 近场 ≥1.5；RCS 单调 | 测量侧一并修 |

预期（诚实版）：G1+G2+G3 之后，近场阶 1.5–1.8（界面场一阶是公认上限，Meep 文档）；若 RCS 与远场量
才是二阶验收对象。到不了 2.0 不一定是错，但"物理误差 vs 真空地板"两栏必须分开报告。

---

## 5. 附录：关键行号

| 主题 | 位置 |
|---|---|
| `ComputeDualGeometry`（scheme/A_stab 签名） | `EmbeddedBoundaryInit.cpp:746-770` |
| 方案 B 的 `Lx/Lz/A_y` | 同文件 `:798-812` |
| `A_stab` 下限 | 同文件 `:825-831` |
| level set 暴露 | `WarpX.H:77`；`AMReX_EBDataCollection.cpp:31-33` |
| parser（局部变量） | `WarpXInitEB.cpp:88-92` |
| one-way / eight-ways XZ 循环 | `WarpXFaceExtensions.cpp:671`、`:831` |
| `GetNeigh/SetNeigh/GetNeighPtr`（dim=0 为 3D 写法） | 同文件 `:277-408` |
| `ComputeSStab` 的 XZ 分支（已就绪） | 同文件 `:96-133` |
| TE 面 BCK-only 特判 | `EmbeddedBoundaryInit.cpp:465-482` |
| `EvolveB` 借用索引（idim=0 无 XZ 分支） | `EvolveB.cpp:295-320`、`:330-350` |
| BCK 兜底 XZ 禁用 | `WarpXFaceExtensions.cpp:612/628` |
| 开关/断言 | `WarpX.cpp:582-596`（上限 2，实现 4 后放开） |

*本文是实施方案，不代替 `C_SIDE_FIX_PLAN.md` 的 F1–F3（已落地）与文档修正；G1 起点即当前 `master`（含 `a55b288` 回退）。*
</content>
