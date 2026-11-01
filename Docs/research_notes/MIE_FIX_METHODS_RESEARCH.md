# 米氏差异的修复方法调研：具体理论 + 操作方法（论文/教科书/开源/本项目源码）

> 日期：2026-10-31 · 配套：`MIE_ERROR_SOURCES_RESEARCH.md`（差异来源拆解）
> 本文对上一份调研列出的六类来源，逐一给出**理论（公式/文献）+ 操作方法（本项目落地步骤/参数/验收）**。
> 结论：真正需要动手的是四件——**① 掩膜改 DM 阈值冻结（理论：Dey–Mittra；操作：`A_y<=A_stab`）；
> ② 精确对偶长度 + B 侧同源（理论：Osher–Fedkiw level set；操作：parser 求交）；③ 近远场变换算 RCS
> （理论：等效原理/Meep `get_farfield`；操作：闭合线 + DFT + Hankel 核）；④ 源/PML 诊断（理论：Oskooi 2008）；
> ⑤ 交错采样一致化（理论：Yee 1966；操作：对理论施加同一平均算子）。**

---

## 0. 总览

| 来源 | 理论出处 | 操作方法（本项目） | 验收数字 |
|---|---|---|---|
| A/F 掩膜/共形/稳定化 | Dey–Mittra 1997；Taflove Ch.8；Benkler 2006 | `A_y <= A_stab` 冻结（DM），激活 `dual_A_y` 通道；标定 `A_stab` | m=0 回 ≥0.85；TE 散射 L2 37.8%→<25% |
| A 几何一阶 | Osher–Fedkiw 2003；AMReX EB2 | parser 精确求交 `Lx/Lz` + B 侧同源积分 | 纯几何阶 ≥1.8；结构残差 24%→<12% |
| B 色散/各向异性 | Taflove Ch.4 | α 分解报告；分辨率扫描；方位平均（仅诊断） | 相位 9–16°→<5°（640²） |
| C 交错采样 | Yee 1966；WarpX `FieldIO` | 对理论施加**同一格心平均算子** | r∈[1,1.2] 误差降一半 |
| D 测量/RCS | Taflove Ch.8（NTFF）；Meep `get_farfield`；Oskooi 2008 | 闭合线 NTFF 直接算 RCS；真空地板诊断 | RCS 与 Mie 差 ≤5%；地板 <5% |
| E Mie 约定 | Bohren–Huffman 1983（级数）；本项目 `mie_cylinder.py` | PEC 内遮罩；scattered 用 `scattered()`；单位核对 | 无（约定问题） |

---

## 1. 掩膜/稳定化：Dey–Mittra 局部共形 + 阈值冻结

### 1.1 理论

- **Dey–Mittra 1997**：对切胞边，把更新系数用"**开放部分面积 / 格面积**"修正，并在轮廓上只取开放段：
  `E_z^{n+1} = E_z^n + (Δt/(ε0 A_open)) ∮_open H·dl`。当 `A_open` 很小时用阈值丢弃（`f < α` → 置零）。
- **Taflove Ch.8.3（"Conformal FDTD"）**：把 DM 系数写成 `T = A_open/A_cell` 的乘积形式，并给出
  "twisted path" 修正；关键结论：**界面场仍一阶**，但切向场/散射截面恢复二阶趋势。
- **Benkler 2006（BCK）**：另一种稳定化——分母抬到 `S_stab`；`S_stab` 由局部最大边长确定。
- **操作含义**：冻结（DM）与抬升（BCK）是两种完整方法；F-3 的教训是"抬到 10 退化阶梯"，
  所以本项目选 **DM 冻结**：`A_y <= A_stab → E=0`。

### 1.2 操作方法（本项目）

1. `MarkUpdateECellsECT` 已支持 `dual_A_y` 可选参数 + `A_stab`；把调用从
   `(m_eb_update_E_fp[lev], edge_lengths_lev)` 改为传 `dual_A_y_lev` 与
   `ect_conformal_E_Astab * dx * dz`（`WarpXInitData.cpp:1679` 处）；
2. `A_stab` 标定：扫描 `A_stab ∈ {0.05,0.1,0.2,0.5,1.0}`（×dx·dz），取"m=0 恢复 + 无热斑 + 环带 L2 最小"；
   预期 **0.1–0.2** 比当前 0.5 更接近 DM 原意（F-3 的 0.5 是为压热斑而设）；
3. 验收：三行判据（全网格 max|E|、PEC 内 max|E|、环带 L2）+ m=0 幅度（目标 ≥0.85）。

### 1.3 文献

Dey & Mittra 1997（DOI 10.1109/7260.622821）；Taflove & Hagness 2005 Ch.8；Benkler et al. 2006（DOI 10.1109/TAP.2006.875494）。

---

## 2. 几何：精确对偶长度 + B/E 同源

### 2.1 理论

- **Osher–Fedkiw 2003**：隐式界面 `φ=0`；沿线段求开放长度 = 对 `[φ≤0]` 的积分；线性插值求交点对线性 φ 精确。
- **AMReX EB2**：每 cell 用节点 φ 构造平面切割（`EB2_2D_C`）；B 侧 `face_areas`/`edge_lengths` 来自同一次构造。
- **同源要求（本项目 E7）**：E 侧若用解析 φ 精确求交，B 侧仍用 AMReX 平面切割 → 两侧几何差 O(h²)（近切点放大为一阶）；
  要真正二阶，必须 **B 侧也换成同一解析 φ 的积分**（face_areas、edge_lengths 重算）。

### 2.2 操作方法（本项目）

1. `ComputeDualGeometry` 的 `scheme=4` 已支持 level-set 线性求交（一阶）；升级为 **parser 精确求交**：
   在 `ComputeDualGeometry` 内用 `eb_implicit_function` 的 parser（`WarpXInitEB.cpp:90` 的同一字符串）逐点求值，
   沿对偶线二分/牛顿求交点（8–16 段足够），得到"近似机器零"的几何误差；
2. **B 侧同源**：在 `ComputeEdgeLengths`/`ComputeFaceAreaFractions` 之后，用同一 parser 重算 `edge_lengths`
   与 `face_areas`（或至少 `edge_lengths`），替换 AMReX 的平面切割值；
3. 验收：纯几何脚本（解析圆）→ `Lx/Lz` 阶 ≥1.8、`A_y` 阶 ≥1.8；L2a 逐位；米氏结构残差 24%→<12%。

### 2.3 文献

Osher & Fedkiw, *Level Set Methods and Dynamic Implicit Surfaces*, 2003；AMReX EB2（`Src/EB/AMReX_EB2_2D_C.cpp`）；
Yu & Mittra 2001（共形介质，同源思想的先例）。

---

## 3. 数值色散/各向异性：定量与报告口径

### 3.1 理论

- **Taflove Ch.4**：Yee 色散 `(1/(cΔt))² sin²(ωΔt/2) = Σ (1/Δx_i²) sin²(k_i Δx_i/2)`；
  沿轴向相位误差 `Δφ ≈ (kΔx)³/24` 每格；**对角方向最差** → 4 重对称伪模（m=4）。
- 量级：20 ppw、~100 格 → ~7.5°；实测 9–16°（含界面与各向异性）。

### 3.2 操作方法

- **α 分解报告**：每次对比给出 `α = Σ E_sim·conj(E_th)/Σ|E_th|²` 的幅度与相位，以及扣 α 后的结构残差；
  把"色散相位"与"几何结构"分开（本次：c0 α=1.008∠9.1°、残差 11%；s2 α=0.977∠16.4°、残差 24%）；
- **分辨率扫描**：40/80/160/640 同口径，看相位是否按 ~h² 收敛；
- **各向异性**：方位模分解（m=4 应随 h² 下降；若不动则是几何/掩膜问题，不是色散）。

---

## 4. 交错采样：对理论施加同一平均算子

### 4.1 理论

- **Yee 1966**：E/B 半格交错；WarpX plotfile 输出是 **cell-centered 平均**
  （`FieldIO.cpp:AverageAndPackVectorField/AverageAndPackScalarField`）：`Ey` 格心值 = 4 个节点值的平均。
- 若理论在格心**点采样**而仿真是 4 点平均，则比较里含平均算子的 O(h²)（光滑区）/O(h)（界面）偏差。

### 4.2 操作方法

- 在 `analyze_mie_convergence.py` / 绘图脚本里，对理论场**做同样的 4 点平均**：
  `Ey_th_cc(i,j) = ¼[Ey(x_i,z_j)+Ey(x_{i+1},z_j)+Ey(x_i,z_{j+1})+Ey(x_{i+1},z_{j+1})]`；
  或直接从 plotfile 读取"节点场"（若 WarpX 支持）避免平均；
- 验收：界面层（r∈[1,1.2]）误差下降 ~一半；不影响远场。

---

## 5. 测量：近远场变换（NTFF）算 RCS，替代"减解析入射"

### 5.1 理论

- **等效原理 / 近远场变换**（Taflove Ch.8.6；Meep `get_farfield`）：在包围散射体的闭合线 `C` 上
  取散射场的 DFT 分量，远场由 Hankel/Green 核积分：
  `S(φ) = ∮_C [ jωμ (n×H_sc)·ŝ + (n×E_sc)·(ŝ×...) ] e^{jk ŝ·r'} dl'`（2D 形式用 H_0^{(1)} 核），
  再 `σ(φ) = (2π/k)|S(φ)|²`（与本项目 `mie_cylinder.rcs` 同口径）。
- **关键**：NTFF 只需要 `C` 上的**散射场**；散射场由"总场 − 入射"得到，入射可以取**数值真空入射**（本项目已有），
  不需要解析 `exp(iKZ)`——这就把"RCS 减解析入射"的误差去掉了。
- **TFSF**（Taflove Ch.5）是更干净的方案：入射只在一层内加减，内部只有散射场；WarpX 无 TFSF，
  用"真空减法"等价替代。

### 5.2 操作方法（本项目）

1. 在 `analyze_mie_convergence.py` 里实现 2D NTFF：取半径 `r_c`（约 1.5R）的闭合线（或矩形），
   用 8 帧 DFT 得到 `E_sc`/`B_sc` 的复振幅（两分量，注意 By/Ey 的 pol 分支）；
2. 入射用**同分辨率真空跑**的 DFT（`A_vac`），`E_sc = A_mie − A_vac`；
3. 2D 远场积分用 `H_0^{(1)}(k|r_c − r'|)`，按角向 DFT 核展开或直接数值积分；
4. 验收：RCS 与 Mie `rcs()` 的 L2 差 ≤5%（当前减解析入射是 19%）；并给"物理误差 vs 地板"两列。

### 5.3 文献

Taflove Ch.8.6；Meep `src/dft.cpp`（`get_farfield`，2D 核实现）；Bohren & Huffman 1983（Mie 级数，RCS 对照）。

---

## 6. 源/PML/地板：诊断与改造

### 6.1 理论

- **Oskooi et al. 2008**：PML 对斜入射、倏逝波、低频有反射；有限厚度形成驻波 → 分辨率无关的残留。
- **Taflove Ch.7**：ABC/PML 的反射系数；源的方向性（CW 激光源有限孔径 → 边缘衍射）。
- 实测：真空地板 14.2%@640² ≈14.7%@160²（分辨率无关）→ 源/PML 伪影。

### 6.2 操作方法

1. **诊断**：分别跑"真空 + 无 PML（增大域）"、"真空 + PML 厚度 ×2"、"真空 + 源孔径变化"，
   看 14% 是否下降；若"无 PML 大域"下降 → PML 反射；若"孔径变宽"下降 → 源衍射；
2. **改造**：加厚 PML（2–4 层 → 8–16 层）、源改为更宽的平面波（或 TFSF，若后续实现）；
   验收：地板 <5%（160²）；
3. 报告：所有 C 侧算例附"真空地板"一列。

### 6.3 文献

Oskooi et al. 2008（DOI 10.1016/j.jcp.2008.01.045）；Taflove Ch.7；Meep PML 文档。

---

## 7. Mie 约定（避免"假差异"）

- **Bohren–Huffman 1983**：Mie 级数的标准形式；本项目 `mie_cylinder.py` 用 2D 圆柱（per-unit-length），
  截断 `n_max = ka + 4ka^{1/3} + 8`；
- **总场 vs 散射场**：`field()` 含入射，`scattered()` 不含；PEC 内 `scattered()` 是解析延拓（不物理，图内遮罩）；
- **单位**：E0 归一化；By 与 Ey 的 B/E 关系（c=1）；RCS 用远场式 `(4/k)|Σ c_n e^{jnφ}|²`；
- 验收：与 `mie_cylinder.py` 自检（`python3 mie_cylinder.py`）对齐。

---

## 8. 落地顺序（与 `MIE_ERROR_SOURCES_RESEARCH.md` 的清单一致）

1. **掩膜 DM 冻结**（§1）→ m=0 恢复、TE 散射 L2 37.8%→<25%；
2. **NTFF RCS**（§5）→ RCS 差 ≤5%；
3. **交错采样一致化**（§4）→ 界面层误差减半；
4. **源/PML 诊断**（§6）→ 地板 <5%；
5. **精确几何 + B 同源**（§2）→ 结构残差 24%→<12%；
6. **色散/各向异性报告口径**（§3）贯穿全程。

---

## 9. 参考文献

1. A. Taflove, S. Hagness, *Computational Electrodynamics: The FDTD Method*, 3rd ed., 2005（Ch.4 色散、Ch.5 TFSF、Ch.7 PML、Ch.8 共形与 NTFF）。
2. K. Yee, IEEE TAP 14(3), 1966（DOI 10.1109/TAP.1966.1138693）。
3. S. Dey, R. Mittra, IEEE MWCL 7(9):273–275, 1997（DOI 10.1109/7260.622821）。
4. S. Benkler, N. Chavannes, N. Kuster, IEEE TAP 54(6):1843–1849, 2006（DOI 10.1109/TAP.2006.875494）。
5. S. Osher, R. Fedkiw, *Level Set Methods and Dynamic Implicit Surfaces*, Springer, 2003。
6. A. Oskooi et al., J. Comput. Phys. 227:5971–5985, 2008（DOI 10.1016/j.jcp.2008.01.045）。
7. C. Bohren, D. Huffman, *Absorption and Scattering of Light by Small Particles*, 1983。
8. Meep：`doc/docs/Subpixel_Smoothing.md`、`src/anisotropic_averaging.cpp`、`src/dft.cpp`（`get_farfield`）。
9. AMReX：`Src/EB/AMReX_EB2_2D_C.cpp`（平面切割）、`AMReX_EBFabFactory.H`（`getLevelSet`）。
10. 本项目：`analyze_mie_convergence.py`、`mie_cylinder.py`、`EvolveE.cpp`、`EmbeddedBoundaryInit.cpp`、`FieldIO.cpp`、`WarpXInitData.cpp`。

*本文是修复方法调研；实测数字与差异拆解见 `MIE_ERROR_SOURCES_RESEARCH.md`。*
</content>
