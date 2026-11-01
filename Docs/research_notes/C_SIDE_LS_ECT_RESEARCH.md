# G1 调研：level-set 精确对偶长度的教科书/论文/开源依据

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 日期：2026-10-31 · 对象：`C_SIDE_LS_ECT_PLAN.md` 的 G1（`Lx/Lz` 的 level-set 积分）
> 目的：确认"沿对偶线求交积分"是标准做法、精度量级、以及 AMReX 侧的数据来源是"同源"的。

---

## 1. 教科书

| 书 | 相关章节/结论 | 对本项目的用法 |
|---|---|---|
| Osher & Fedkiw, *Level Set Methods and Dynamic Implicit Surfaces*, Springer 2003 | 用符号距离函数 `φ` 隐式表示界面；界面位置由 `φ = 0` 给出；**线性插值求交**在网格分辨率下是标准的二阶做法 | `Lx/Lz` 的积分算法（采样 + 符号变化 + 线性插值） |
| Taflove & Hagness, *Computational Electrodynamics*, 3rd ed. 2005 | 共形 FDTD 的几何修正必须用"开放面积/开放边长"，且二者需成对 | 为什么只改 `Lx/Lz` 而不改 `A_y` 是自洽的（`A_y` 已是面积分） |
| LeVeque, *Finite Volume Methods for Hyperbolic Problems*, CUP 2002 | cut-cell 的守恒量是**控制体内的积分**，边界用精确几何截断 | `⅓` 式的一维截断积分（不是分数平均）的合法性 |
| Bossavit, *Computational Electromagnetism*, 1998 | 离散 Hodge star 的精度取决于几何量的积分精度 | 支撑"几何量必须与真几何同阶"的判断 |

**关键点**：Osher–Fedkiw 的"用 `φ` 的线性插值定位界面"正是本实现的核心——`φ` 在 cell 内被 EB2
线性化为平面，因此**沿线的线性插值求交对该平面是精确的**，对真曲面是 `O(h²)`。

## 2. 论文

| 文献 | 结论 | 引用位置 |
|---|---|---|
| Dey & Mittra 1997（IEEE MGWL 7(9):273–275） | 局部共形 FDTD 的原始定义：Faraday 用**开放面积**、Ampere 用**开放长度/面积**，两者必须成对 | G1 的物理依据 |
| Farjadpour et al. 2006（Opt. Lett. 31:2972） | Meep 的 subpixel smoothing：**用 level set 解析计算**每个 cell 的填充分数与法向，而不是采样体积分数 | "level set → 几何分数"的权威先例；说明 `volfrac` 是 `φ` 的**一阶近似** |
| Kottke, Farjadpour, Johnson 2008（PRE 77:036611） | subpixel 的误差分析：直接阶梯化把二阶降为一阶 | 解释为什么当前的 `½(v+v)` 是一阶 |
| Nieter et al. 2009（JCP 228:7902） | VORPAL 的 Dey–Mittra 实现；小胞元使 `dt` 任意缩小 | 提醒 level-set 之后仍需稳定化（G3） |

**结论**：Meep 用的是"level set → 解析填充分数"，正是 G1 要做的；它同时说明**为什么要用 `φ` 而不是
`volfrac`**——`volfrac` 是 `φ` 在 cell 上的积分，用它反推"线长"必然损失一维信息（维度错配）。

## 3. 开源项目（源码级核对）

| 项目 | 事实 | 出处 |
|---|---|---|
| **AMReX** | `EBFArrayBoxFactory::getLevelSet()` 返回 `const MultiFab&`，`IntVect::TheUnitVector()`（**全 nodal**）、1 comp、`nGrow = m_ngrow[0]`；由 `a_level.fillLevelSet(*m_levelset, m_geom)` 填充 | `AMReX_EBFabFactory.H:92`、`AMReX_EBDataCollection.cpp:31-33` |
| **AMReX** | `volfrac` 是 cell-centered、由同一次 `fillEBCellFlag`/几何构造得到 → **与 level set 同源** | 同上 |
| **WarpX** | `WarpX::fieldEBFactory(lev)` 暴露该工厂；`WarpX.H` 通过 `m_ebdc` 持有 | `WarpX.H:77` |
| **Meep** | `geom.py` 用 level set 的解析表达式算 subpixel 平均（不是采样） | `doc/docs/Subpixel_Smoothing.md` |

**同源性结论**：`getLevelSet()` 与 `getVolFrac()`/`getAreaFrac()` 来自**同一次 EB 构造**，
所以用 level set 算 `Lx/Lz` 不会引入"E 侧与 B 侧不同源"的新问题（E7 风险不增加）。

## 4. 精度量级（预期）

| 环节 | 误差 | 说明 |
|---|---|---|
| EB2 的 `φ` 对真曲面的表示 | `O(h²)`（每 cell 线性化） | 这是 AMReX 的固有精度，B 侧也一样 |
| 沿线的线性插值求交 | 对线性 `φ` **精确**；对曲 `φ` 为 `O((h/N)²)` | `N = 8` 时 `(h/8)² ≈ 0.016h²` |
| 采样求和 | `O(h/N)` 的分段常数误差 | 用"符号变化子区间 + 线性插值"消掉 |

**总计 `O(h²)`**，与 `A_y` 同阶。这就是 G1 能把几何阶数从 ~1.0 提到 ~2 的原因。

## 5. 与方案 B 的对比（为什么必须换）

| | 方案 B（现状） | 方案 C（G1） |
|---|---|---|
| 输入 | `volfrac`（2D 面积分数） | `φ`（nodal 标量场） |
| 运算 | 相邻两个 cell 的 `v` 取平均 × h | 沿对偶线求 `φ ≤ 0` 的长度 |
| 维度 | 2D → 1D（**信息损失**） | 1D → 1D（维度匹配） |
| 实测阶 | 1.02 / 0.76 | 预期 ≥ 1.8 |

## 6. 实现要点（从调研直接导出）

1. **符号**：`φ > 0` 为固体（WarpX 的 `eb_implicit_function` 约定），开放是 `φ ≤ 0`；
2. **z 向平均**：`Lx` 在 `z = z_{j+½}`，而 `φ` 在节点 → `φ_c = ½(φ(i,j) + φ(i,j+1))`（z 向线性）；
3. **x 向插值**：采样点 `x = x_i − h/2 + s·h/N`，落在 `[x_{i−1}, x_i]` 或 `[x_i, x_{i+1}]`，
   用该段的两个节点线性插值；
4. **规则域短路**：若 4 个角点 `φ` 全部 `≤ 0`，直接返回 `h`（保证 L2a 逐位）；
5. **完全不覆盖**：全 `> 0` → 返回 0；
6. **不改变 `A_y`**：仍用 `¼h²Σv`（已是二阶）；
7. **仍施加 `A_stab` 下限**（G3 之前保持 0.1），但 `A_y == 0` 保持 0。

## 7. 参考文献

1. S. Osher, R. Fedkiw, *Level Set Methods and Dynamic Implicit Surfaces*, Springer, 2003.
2. A. Taflove, S. C. Hagness, *Computational Electrodynamics*, 3rd ed., Artech House, 2005.
3. R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, CUP, 2002.
4. A. Bossavit, *Computational Electromagnetism*, Academic Press, 1998.
5. S. Dey, R. Mittra, IEEE Microwave Guided Wave Lett. **7**(9), 273–275, 1997.
6. A. Farjadpour et al., *Opt. Lett.* **31**(20), 2972–2974, 2006.
7. A. Kottke, A. Farjadpour, S. G. Johnson, *Phys. Rev. E* **77**, 036611, 2008.
8. C. Nieter et al., *J. Comput. Phys.* **228**(21), 7902–7916, 2009.
9. AMReX: `Src/EB/AMReX_EBFabFactory.H:92`、`Src/EB/AMReX_EBDataCollection.cpp:31-33`。
10. WarpX: `Source/WarpX.H:77`（`fieldEBFactory`）、`Docs/source/glossary.rst`（ECT/BCK）。
