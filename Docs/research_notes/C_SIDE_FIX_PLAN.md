# C 侧修正方案（依据 `C_SIDE_REPORT_CHECK.md`）

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 日期：2026-10-31 · 输入：`C_SIDE_REPORT_CHECK.md`（第三方复核，复跑 L2a/L3/米氏 + 纯几何实验）
> 输出：本文只给**思路与改动清单**，不改代码。
> 涉及文件：`Source/WarpX.H`、`Source/WarpX.cpp`、`Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp`、
> `Source/EmbeddedBoundary/EmbeddedBoundaryInit.{H,cpp}`、`Source/Initialization/WarpXInitData.cpp`、
> `benchmarks/cases/c02_0_pec/analyze_mie_convergence.py`、`HANDOFF.md`、`C_SIDE_REPORT.md`。

---

## 0. 对 CHECK 的逐条处置

CHECK 的复跑结论（数字全部复现、报告可信）我接受，不辩解。它指出的偏差分三类，处置如下：

| CHECK 项 | 判定 | 处置 |
|---|---|---|
| **P1：R5 机理写错**（初值 1.0 → 多层低层 ÷1.0 → **量纲错但不崩**，不是我写的"除零保护 → 冻结"） | **接受，我错了** | 改 `C_SIDE_REPORT.md` §4.1 R5 的描述；代码加多层断言（见 §2.2） |
| P3：`WarpX.H` 注释写 "from the EB level set"，实现是 volfrac | **接受** | 改注释，并把方案编号显式化 |
| P3：掩膜路径 `update==0` 只 `return`、不写 0 | **接受（补强）** | 改成写 0；同时把"covered 强制 0"的语义补齐 |
| P3：`ect_conformal_E` 无下限断言，负数静默关闭 | **接受** | 加 `0 <= x <= 4` 断言 |
| P2：`HANDOFF.md` §3/§4 过期（C 侧仍写"未做"、归因停在"CW 源"） | **接受** | 改 §3 为"多 grid"、§4 删掉"C 侧未做" |
| §5.1：`Lx/Lz` 一阶（RMS 0.097→0.048→0.027，局部可达 0.5h），`A_y` 二阶（1.97） | **接受，且这是本轮最有价值的发现** | 直接决定方案 A/C 的优先级（见 §1.2） |
| §5.1：frozen-but-open 节点 176/356/702（≈5.6 个/单位周长） | **接受** | 掩膜改 `A_y <= 0`（见 §2.3），改后重测该统计 |
| §5.2：真空平面波残差 83%/21%/15%；RCS 减的是解析入射 | **接受** | 测量侧修正（见 §2.4），并把真空残差当误差地板报告 |
| §6 建议"先改 R5 描述" | 接受 | 放在第一步（零风险） |

**唯一需要补充的一条**（CHECK 没有明确区分）：它把"修正版方案 A（面分数）"和"方案 C（level set）"
并列为"才有希望接近 2 阶"的路径。我认为两者的量级不同，必须分开验收：

- **方案 A（面分数，四面平均）**：把"用 2D 面积分数估 1D 长度"换成"用 1D 面长度分数估 1D 线长度"，
  消除了**维度错配**这个主导误差，但面分数与格心线仍差半格 → 对曲边界**仍是一阶**，只是系数更小；
- **方案 C（level set 求交）**：把几何量算成**精确的开放长度/面积**，才可能到二阶。

所以验收判据要分开：方案 A 只要求"Lx/Lz 的 RMS 误差显著下降（例如 0.097→0.05 以下）"，
方案 C 才要求"RMS 阶数 ≥ 1.8"。

---

## 1. 修正思路

### 1.1 三层，按"风险 × 收益"排序，不混在一起

```
L1 文档/注释纠错         零风险，先做，避免错误结论继续传播
L2 语义与防护           小改动，把"静默错误"变成"显式错误"
L3 几何与测量升级        真正决定能否到二阶，需要实验迭代
```

**为什么必须分层**：L3 会改变所有判据数字（L2a 之外的基线全部要重记）。如果 L1/L2 混在 L3 里，
一旦 L3 结果不理想，就分不清是几何改坏了还是原本就有问题。

### 1.2 L3 为什么先改 `Lx/Lz` 而不是先上 level-set

CHECK 的几何实验给了明确排序：

| 量 | 当前实现 | 绝对 RMS 阶 | 性质 |
|---|---|---|---|
| `A_y` | `¼dxdzΣvolfrac`（4 象限） | **1.97** | 已经是二阶，不需要改 |
| `Lx`/`Lz` | `½dx(v+v)`（方案 B） | **1.0** | 维度错配（2D 面积分数 → 1D 长度） |

近切点处 `Lx ≈ 0.5h`（真实值 ≈0），这个"半格"误差直接进入
`Ey` 的通量权重和 `Ex`/`Ez` 的分母，是当前一阶误差的主源。

因此顺序是：**先把 `Lx/Lz` 从"面积分数"换成"面分数"（方案 A，维度正确）→ 测几何 RMS 是否降一个量级
→ 若仍不够再上 level-set（方案 C）**。这样能用半天时间判断"维度错配"是否就是主导项，
避免直接投入 level-set 的实现成本。

### 1.3 掩膜必须同步改，否则几何改了也用不上

CHECK 统计出 frozen-but-open 节点 176/356/702（≈5.6 个/单位周长）。这些节点 `A_y > 0` 却被
`ly == 0` 冻结，几何再好也不参与更新。**几何与掩膜是同一个改动的两半**，必须一起改、一起验。

### 1.4 测量侧要和物理侧分开报告

真空残差 15%（160²）意味着：即使几何完全精确，米氏 L2 也降不到 0.15 以下。所以修正后的报告要写成

```
L2_total ≤ L2_physical + L2_floor(真空残差)
```

否则会把"地板"当成"格式还没到二阶"。RCS 同理：必须减**同分辨率的数值入射**，而不是解析 `exp(iKZ)`。

---

## 2. 改动清单（逐文件、逐位置）

### 2.1 L1：文档与注释（零风险，第一步）

| 文件 | 位置 | 改什么 |
|---|---|---|
| `Docs/research_notes/C_SIDE_REPORT.md` | §4.1 表格 R5 行 | 把"低层几何为 0 → 除零保护 → E 全冻结"改为：**分配初值 1.0，`ComputeDualGeometry` 只在 `maxLevel()` 调用，低层保持 1.0，内核 ÷1.0 → 量纲错、静默错误更新（不冻结、不 abort）** |
| `Source/WarpX.H` | `:198-201` 注释 | `2 = on, with the real dual geometry from the EB level set` → `2 = on, dual geometry assembled from EB volume fractions (scheme B); 3 = face fractions (scheme A); 4 = level set (scheme C, TODO)` |
| `Docs/research_notes/HANDOFF.md` | §3 表格、§4 待办 1 | §3 补一行"多 grid"结论；§4 把"C 侧（未做）"改为"C 侧已实现（tag `c-side-2d-v1`），待做 R1/R2/R3" |
| `CASES_2026-10-31.md` | C16 结论段（§7sexies.6） | 追加一句 R5 的正确机理（与 R5 一致） |

### 2.2 L2：语义与防护（小改动，一次编译可验）

| 文件 | 位置 | 改什么 |
|---|---|---|
| `Source/WarpX.cpp` | `:582-596` | 断言从 `ect_conformal_E == 1 \|\| == 2` 改为 `0 <= ect_conformal_E <= 4`；负数/超范围立即 abort |
| `Source/WarpX.cpp` | `:2714-2727` 分配段 | 在 `if (ect_conformal_E > 0)` 内加 `WARPX_ALWAYS_ASSERT_WITH_MESSAGE(maxLevel() == 0, "algo.ect_conformal_E does not support mesh refinement yet")` |
| `Source/FieldSolver/FiniteDifferenceSolver/EvolveE.cpp` | `:290` / `:301` / `:313` | 掩膜分支从 `{ return; }` 改为 `{ Ex(i,j,k) = 0; return; }`（三分量同理），把"covered 强制 0"落实 |
| `Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp` | `:798-812` | 新增 `scheme` 参数（1=B 常数、2=B、3=A、4=C），`Lx/Lz` 按方案分支；`A_y` 保持不变（已是二阶） |

**方案 A 的实现要点**（`scheme == 3`）：改用 `getAreaFrac()` 的两个面（AMReX 的 lo 面），四面平均：

```cpp
auto const& ax = area_frac[0]->const_array(mfi);   // x-lo 面：沿 z 的开放长度分数
auto const& az = area_frac[1]->const_array(mfi);   // z-lo 面：沿 x 的开放长度分数

Lx(i,j,k) = 0.25_rt*dx*( az(i-1,j,  k) + az(i-1,j+1,k)
                       + az(i,  j,  k) + az(i,  j+1,k) );
Lz(i,j,k) = 0.25_rt*dz*( ax(i,  j-1,k) + ax(i,  j,  k)
                       + ax(i+1,j-1,k) + ax(i+1,j,  k) );
```

索引必须用 S1a' 的单 cell 真值表再核一次（`getAreaFrac` 是 lo 面，`i-1`/`i+1` 的归属容易差一格）。

### 2.3 L3-a：掩膜改成 `A_y <= 0`

| 文件 | 位置 | 改什么 |
|---|---|---|
| `Source/Initialization/WarpXInitData.cpp` | `:1636-1655` | 调整调用顺序：先 `ComputeDualGeometry`，再 `MarkUpdateECellsECT`，并把 `dual_A_y` 传进去 |
| `Source/EmbeddedBoundary/EmbeddedBoundaryInit.H` | `:75` | `MarkUpdateECellsECT` 增加 `amrex::MultiFab const* dual_A_y` 形参（默认 `nullptr`，保持旧行为） |
| `Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp` | `:233-292` | XZ 的 `Ey` 分支：若 `dual_A_y` 非空，判据改为 `A_y(i,j) <= 0`；否则沿用 `ly` |

验收：改后打印"被冻结但 `A_y > 0`"的节点数，应为 0（当前 176/356/702）。

### 2.4 L3-b：测量修正（不改物理，只改判据）

| 文件 | 位置 | 改什么 |
|---|---|---|
| `benchmarks/cases/c02_0_pec/analyze_mie_convergence.py` | `:118-140`（`Asc = A_al - exp(iKZ)`） | 改为减**同分辨率真空跑的数值入射**：新增 `--vacuum <dir>` 参数，读 `vac/ey_<n>` 的谐波幅值 `A_vac`，用 `Asc = A_al - A_vac`（两者已用同一全局相位对齐） |
| 同上 | 输出段 | 把真空残差（`‖A_vac − exp(iKZ)‖/‖exp(iKZ)‖`）作为 `floor` 一并打印，并给出 `L2_phys = sqrt(max(L2² − floor², 0))` |
| `make_mie_cases.py` | `NCELLS` | 增加 320 档（至少 4 档再拟合阶数） |

### 2.5 L3-c：方案 C（level set），最后做

| 文件 | 位置 | 改什么 |
|---|---|---|
| `Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp` | 新函数 `ComputeDualGeometryFromLevelSet` | 取 `eb_fact` 的 level set（`EB2::Level::getLevelSet()`），沿格心线采样 8–16 点求 `φ <= 0` 的段长（**注意本项目约定：正值为固体**），积分得 `Lx/Lz`；`A_y` 仍用 4 象限精确面积 |
| 同上 | `scheme == 4` 分支 | 调用上述函数 |

验收：几何 RMS 阶数 ≥ 1.8；米氏近场阶数在扣掉真空地板后 ≥ 1.7。

---

## 3. 执行顺序与验收（每步单独提交）

| 步 | 内容 | 验收判据 | 风险 | 工时 |
|---|---|---|---|---|
| **F1** | §2.1 文档/注释 | 无数字变化；`git grep "level set" Source/WarpX.H` 不再出现在 `=2` 的描述里 | 无 | 0.2 h |
| **F2** | §2.2 断言 + 掩膜写 0 + 多层 abort | 编译通过；L2a 仍 PASS（rel ≤ 1e-14）；`ect_conformal_E = -1` 立即 abort | 低 | 0.5 h |
| **F3** | §2.3 掩膜改 `A_y <= 0` | 冻结但开放节点数 = 0；L2a 仍 PASS；L3 45° 的 `\|r+1\|` 记录新值 | 中（顺序调整容易漏 `FillBoundary`） | 2 h |
| **F4** | §2.2 新增方案 A（`scheme=3`），作为默认 | 纯几何 RMS：`Lx/Lz` 从 0.097/0.048/0.027 降到 ≤0.05/0.025/0.012（约一半）；`A_y` 不变（1.97 阶） | 中（lo 面索引） | 4 h |
| **F5** | 重跑米氏 4 档 × 2 开关，记录 `L2`、真空地板、`L2_phys` | 报告"物理误差 vs 地板"两栏；阶数不再由地板主导 | 低（但要算力） | 2 h |
| **F6** | §2.4 RCS 减数值入射 | RCS 阶数在 4 档上单调；若仍 1.7，说明是角度平均的固有平滑 | 低 | 2 h |
| **F7** | §2.5 level-set（方案 C） | 几何 RMS 阶 ≥ 1.8；`Lx/Lz` 在近切点不再有 0.5h 误差 | 高（level set 采样 + 符号约定） | 8 h |

**总工时约 19 h（2.5 个工作日）**，其中 F1–F3 是"纠错 + 让已有几何真正生效"，F4–F7 才是"冲二阶"。

---

## 4. 如果 F4 之后阶数仍不到 1.7，怎么办

按可能性排序，逐条排除（不要同时改两个）：

1. **掩膜没真正放开**：打印 `A_y<=0` 与 `update==0` 的节点集合，确认后者是前者的子集；
2. **`Lx/Lz` 的 lo 面索引差一格**：用单 cell 真值表（S1a'）复核 `a_x(i+1,j)` 的归属；
3. **B 侧几何与 E 侧不同源**（RESEARCH2 §3.2）：`face_areas` 来自 `area_frac`，而 E 侧用四面平均，
   两者对同一 cell 的"开放"描述不一致 → 需要把 B 侧也换成同样的平均（工作量 +1 天）；
4. **界面场本身至多一阶**（Meep 文档的 limitation 4）：那近场判据应移动到更外层（如 `3R–5R`），
   或在 RCS 上做验收，而不是死磕近场；
5. **真空地板仍主导**：把 320 档也做真空对照，确认地板随 `dx` 下降（应该 ~1 阶，因为色散是二阶但
   拟合窗固定会引入一阶项）。

---

## 5. 不做/暂缓（附理由）

| 项 | 决定 | 理由 |
|---|---|---|
| 对偶稳定化 `A_stab` | 暂缓到 F4 之后 | 当前算例（单 grid）没有小胞元爆炸；先拿到干净的几何基线，否则稳定化会把精度信号糊掉 |
| J 项对偶积分 | 暂缓 | 当前 `J = 0`；加粒子时再改，且需要 Esirkepov 型沉积配合（M2a-2） |
| PML 区排除 | 暂缓但**加断言** | 现有算例 EB 都在内部；加"EB 与 PML box 相交则 abort"比实现区排除便宜 |
| F 场（divE 清洁） | 暂缓但**加断言** | `Ffield != nullptr` 且开关打开时 abort，避免静默丢项 |
| 3D / AMR | 排在 F7 之后 | 依赖 2D 的几何与掩膜方案定型；AMR 还要 per-level 初始化 |

---

## 6. 一句话总结

CHECK 让我把"残留问题"里最含糊的一条（R5 机理）纠正了，并且用几何实验证明：
**当前的 1.06 阶不是"共形方法不行"，而是 `Lx/Lz` 用 2D 面积分数估 1D 长度（维度错配，局部半格）
+ 掩膜冻结 5.6 节点/周长 + 真空地板 15% 三者叠加**。

所以修正顺序是：**先纠错（F1）→ 再让已有几何生效（F2/F3）→ 再换维度正确的面分数（F4）→
最后才上 level-set（F7）**，每一步都有独立的、可证伪的验收判据，不把"能不能到二阶"押在单次大改上。

---

## 7. 方案修订（执行前收到的 6 点，2026-10-31）

以下 6 条全部采纳，并替代上文中与之冲突的表述。

| # | 问题 | 采纳后的做法 |
|---|---|---|
| **M1** | F2 把断言放宽到 `0..4`，但 scheme 3/4 尚未实现；`constant_geometry` 是 bool，非 1 的值会被当成"真实几何 B"，**静默回退** | F2 阶段断言**只允许 `0..2`**；同时把 `constant_geometry` 改成显式 `scheme` 参数（1=常数、2=B），`3/4` 立即 `ABORT("not implemented")`。放宽到 3 的那一步与 F4 的代码同一个提交 |
| **M2** | F4 若把方案 A 定为默认，会一次性改掉全部基线 | F4 **保留默认 2**，用显式 `=3` 验证；F5/F6 通过后再单独一个提交改默认并统一重记基线。**L2a 回归加入 F4/F7 的验收**（此前只有 F2 有） |
| **M3** | `Asc = A_al − A_vac` 的相位：`A_al` 已乘过"对 Mie 的全局相位"，`A_vac` 若不再乘同一因子，RCS 会多一个整体旋转；`sqrt(L2²−floor²)` 是正交近似 | `A_vac` 乘**同一个** `exp(i·phase)`；报告里**直接给 `L2` 与 `floor` 两列**，`L2_phys` 只作辅助不宣称"物理阶数" |
| **M4** | "掩膜写 0"只改了共形内核，`EvolveE.cpp` 的 Yee 内核（`:175/186/197`）仍是 `return` | 只改共形内核（`:292/301/310`），**不动 Yee**（改了等于改上游行为）；在代码注释与本文里注明"两套内核的掩膜语义不同" |
| **M5** | 方案 A 的 `az(i−1,j+1)` / `ax(i+1,j−1)` 在最外层 grow 会超出 `area_frac` 的 ghost | 写入 box 限制为 **`valid + 1` 层**（不是 `nGrowVect()`），并对 `i−1/j+1` 做 clamp，与 `ComputeEdgeLengths` 的 `lx_lo/lx_hi` 做法一致 |
| **M6** | §2.1 引用 `CASES_2026-10-31.md §7sexies.6` 不够明确 | 统一写成 `CASES_2026-10-31.md` 的 **C16 结论段（§7sexies.6）** |

对应的执行顺序调整为：

```
F1  文档/注释纠错（含本文 §7）
F2  断言 0..2、多层 abort、scheme 参数化、共形内核掩膜写 0   → 验：-1/3 必须 abort；L2a PASS
F3  掩膜改 A_y<=0                                            → 验：冻结但开放节点数 = 0；L2a PASS
F4  方案 A（scheme=3），默认仍为 2                            → 验：几何 RMS 减半；L2a(=3) 仍 PASS
F5  米氏 4 档 + 真空地板                                     → 验：两列报告
F6  RCS 减数值入射（同一 phase 因子）                         → 验：阶数单调
F7  level-set（scheme=4），之后才考虑改默认                   → 验：RMS 阶 >= 1.8；L2a(=4) 仍 PASS
```
