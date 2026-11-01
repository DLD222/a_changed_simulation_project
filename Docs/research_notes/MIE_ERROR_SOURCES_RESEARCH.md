# 米氏对比中"仿真 vs 理论"差异的来源调研（论文/教科书/开源/本项目源码）

> 日期：2026-10-31 · 对象：2D TE/TM 平面波照射 PEC 圆柱的 Mie 对照（`benchmarks/cases/c02_0_pec`）
> 实测基线（160²，散射场 L2 相对误差）：**TE 37.8%、TM 19.7%**；c0（掩膜 Yee E）近似 32.7%（总场）。
> 本文回答一个问题：**这些差异分别来自哪里**——按"几何/界面、色散、采样后处理、测量约定、级数约定、本项目特有"六类，
> 每类给出文献依据、本项目源码位置、量化估计和实测签名。

---

## 0. 结论摘要

| 来源 | 文献/开源依据 | 本项目源码位置 | 预期阶/量级 | 实测签名（本次） |
|---|---|---|---|---|
| A 几何/界面表示（阶梯→共形、对偶一阶、掩膜冻结） | Dey–Mittra 1997；Yu–Mittra 2001；Benkler 2006；Xiao–Liu 2004/08；Taflove Ch.8 | `EvolveE.cpp`（ECT 内核、掩膜）、`EmbeddedBoundaryInit.cpp`（A_stab、Lx/Lz） | 界面场一阶；几何 1.0 阶 | m=0 幅度低（s2 0.65 vs c0 0.86）；结构残差 24%（s2）vs 11%（c0） |
| B 数值色散 + 网格各向异性 | Taflove & Hagness 2005 Ch.4；Meep subpixel 文档 | 全局（Yee 交错 + 中心差分） | 相位 ~7–20°；各向异性 ~4 重对称 | 相位 9–16°；m=3/m=4 异常 |
| C 交错采样与后处理平均 | Yee 1966；WarpX `FieldIO` 平均 | `FieldIO.cpp:AverageAndPack*` | 光滑区 O(h²)，界面 O(h) | r<1.1 附近误差抬升（TE 75%@1.05） |
| D 测量/基准约定（DFT、对齐、RCS 减解析入射、真空地板 14%、PML/源） | Oskooi 2008 PML 综述；Taflove Ch.7 | `analyze_mie_convergence.py` | 地板分辨率无关 | 真空 14.2%（640²）≈14.7%（160²） |
| E Mie 级数/单位约定 | 本项目 `mie_cylinder.py` 文档串 | `mie_cylinder.py` | 截断误差极小 | PEC 内为解析延拓（图内不该画） |
| F 本项目特有（A_stab、掩膜、BCK、借用、无 F 项） | 见 §6 | `EvolveE.cpp`、`a55b288`/`6b76d4a` | 主要 A/F 项 | 与 c0 的差 = 2× 结构残差 |

**一句话**：TE 的 37.8% 里，最大一块是 **F（本项目 E 侧掩膜/稳定化削弱表面电流）** 与 **A（几何一阶）**；
其次 **B（相位/各向异性）**；**D（14% 地板）** 决定 640² 之后的下限；TM 好得多是因为 Ex/Ez 的掩膜判据（lx/lz）
比 Ey 的 `ly` 判据温和。

---

## 1. 差异的定量分解（先把"差异"拆开）

### 1.1 全局复因子 vs 结构残差（散射场，r∈[1.1,2.5]）

| 配置 | 全局因子 α（幅度∠相位） | 扣 α 后结构残差 / 场 |
|---|---|---|
| c0（掩膜 Yee E） | 1.008 ∠ +9.1° | **11.2%** |
| s2（C 侧 volfrac） | 0.977 ∠ +16.4° | **24.1%** |

- **c0 的误差主要是"全局相位 9°"**（色散 + 源），结构误差只有 11%；
- **s2 的相位偏移加倍（16°），结构残差翻倍（24%）**——这就是 C 侧"目前还不如 c0"的量化解释。

### 1.2 径向分布：误差不是贴壁热斑

| r | TE 误差 | TM 误差 |
|---|---|---|
| 1.05 | 75.3% | 20.8% |
| 1.10 | 32.4% | 19.6% |
| 1.40 | 35.6% | 20.2% |
| 2.00 | 37.8% | 21.0% |
| 3.00 | 40.1% | 22.2% |

r≥1.1 之后几乎**平坦**——说明差异是**全局模态误差**（幅度/相位/各向异性），不是贴壁局部伪影。

### 1.3 方位模分解（TE 散射场 |c_m|，r=1.5）

| m | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| c0 sim/theory | 0.86 | 0.94 | 1.02 | 1.13 | 1.10 | 1.01 | 0.94 |
| s2 sim/theory | **0.65** | 1.18 | 0.89 | **2.16** | 1.13 | 0.93 | 0.88 |
| s2 相位(°) | 25 | 15 | 23 | 33 | 26 | 10 | −1 |

- **m=0（单极/前向阴影）被压低**：s2 只剩 65%（c0 86%）→ 表面电流被削弱（★ F 类，见 §6）；
- **m=3 异常放大**（理论该模极小，sim 有 2–9×）→ 方形网格的**各向异性伪模**（★ B 类）；
- **m=4 ≈1.1**（c0/s2 都有）→ 4 重对称的网格各向异性（★ B 类）；
- 相位整体 +10–30° → 色散 + 掩膜相位（★ A/B 类）。

---

## 2. 来源 A：几何/界面表示（阶梯、共形、对偶几何、掩膜）

### 2.1 文献

- **Dey & Mittra 1997**（Locally conformal FDTD）：阶梯面把曲面 PEC 变成"面内偏移"，误差一阶；
  共形技术用局部修正的场更新把边界"贴"到真曲面，但**只对切向场二阶、界面场仍一阶**。
- **Yu & Mittra 2001**（共形介质曲面）：共形面积/长度的等效系数，界面上二阶、内部一阶的混合。
- **Benkler, Chavannes & Kuster 2006**（BCK）：小面积面的稳定化（抬到 S_stab）——稳定但引入局部偏差。
- **Xiao & Liu 2004/2008**（ECT）：扩展/借用——把邻面电动势借给不稳定面；同样只解决稳定性与一阶几何。
- **Taflove & Hagness 2005**（Ch.8 "Conformal FDTD"）：综述阶梯→共形的误差阶。

### 2.2 本项目源码与实测对应

| 机制 | 源码 | 对差异的贡献 |
|---|---|---|
| B 侧 ECT（面面积/借用/BCK） | `EvolveB.cpp` + `WarpXFaceExtensions.cpp` | 曲面电流的"生成面"精度；m=4 各向异性的一部分 |
| E 侧对偶几何 | `EvolveECartesianECT`（`Ex/Lz`、`Ez/Lx`、`Ey/A_y`） | 界面场一阶 → m=0 幅度与相位 |
| `Lx/Lz` 构造（scheme 2 volfrac 或 4 level-set 线性求交） | `EmbeddedBoundaryInit.cpp` | 纯几何实测 **~1.0 阶**（0.029/0.014/0.008 m）→ 结构残差 |
| `A_stab` 下限（默认 0.5·dx·dz） | `a55b288`、`6b76d4a` | 小对偶面被抬升 → 局部欠松弛 |
| **`ly` 掩膜冻结切面节点** | `MarkUpdateECellsECT`（调用处仍传 `ly`） | **表面电流被冻结节点截断** → m=0 掉 35%（与 c0 的 14% 对比即此项） |
| `EvolveE.cpp` 中 `update_*==0 → return` | 内核 | 同上（冻结不参与更新） |

**结论**：A 类（几何一阶 + 掩膜冻结）是 TE 误差的第一大项；其中"掩膜"又是 **C 侧特有**（c0 没有这层冻结）。

---

## 3. 来源 B：数值色散与网格各向异性

### 3.1 文献

- **Taflove & Hagness 2005, Ch.4**：Yee 格式的数值色散关系；二阶中心差分的相位误差
  `Δφ ≈ (kΔx)³/24` 每格（沿轴向），且**依赖传播方向**（对角方向色散最强）→ 4 重对称伪模（m=4）。
- **Taflove Ch.4 的"网格各向异性"**：方网格对圆柱散射会引入 m=4 误差；细网格按 O(h²) 收敛。
- **Meep 文档（Subpixel_Smoothing.md）**：明说"sharp corners/edges 产生场奇点、降低误差阶"，
  且界面场至多一阶——与我们的 m=3/m=4 伪模一致。

### 3.2 量级估计

- 160²：λ/Δx = 20 ppw；每格相位误差 ≈ (kΔx)³/24 = (0.314)³/24 ≈ 1.3e-3 rad；
  跨 100 格（约 5λ）≈ 0.13 rad ≈ **7.5°**；加上各向异性与界面，实测 9–16° 自洽。
- 640²：相位误差降 64× → 相位应趋零；但**总差只剩 ~18–22%**，说明此时其余项与**测量地板（14%）**主导。

### 3.3 实测签名

- 低阶模相位 +10–30°、随 r 平缓积累（典型色散特征）；
- m=4 ≈1.1（c0/s2 都有）；m=3 在 r=2.0 处被放大到 2–9×（理论该模仅 5.8e-3，属"小模被伪模淹没"）。

---

## 4. 来源 C：交错网格采样与后处理平均

### 4.1 文献

- **Yee 1966**：E/B 半格交错；理论公式是连续解，差分场在交错点上。
- 本项目 **plotfile 输出是 cell-centered 平均**：`Source/Diagnostics/FieldIO.cpp` 的
  `AverageAndPackVectorField/AverageAndPackScalarField` 把交错场平均到格心——`Ey`（节点量）
  在格心 = 4 个节点的算术平均。

### 4.2 影响

- 光滑区：平均是 O(h²) 算子，误差小；
- 界面附近（r≈1）：平均跨过 PEC 边界 → O(h) 偏差；这正是 r=1.05 处 TE 误差 75% 的一部分；
- 与 C 侧无关（c0 同款），不是主要项，但应在界面层报告时注明。

---

## 5. 来源 D：测量/基准约定（本项目脚本）

### 5.1 现状（源码）

| 项 | 位置 | 影响 |
|---|---|---|
| 8 帧 DFT | `analyze_mie_convergence.py`（`load_frames`/`harmonic_amplitude`） | 窗口含瞬态/两频泄漏 → 幅度/相位偏差；实测 DFT/LSQ 差 1.5×（40²） |
| 全局相位对齐（总场 vs Mie 总场） | 同上 | 对齐残差直接进入相位误差 |
| RCS 减**解析** `exp(iKZ)` | 同上（RCS 分支） | 数值入射 ≠ 解析入射 → RCS 偏差；本次已用真空数值入射做散射场 |
| 近场环带 1.5R–2.5R | 同上 | 排除界面层；但 15% 地板仍在 |
| 真空地板 | `/tmp/mie640/vac` 实测 | **14.2%@640² ≈14.7%@160²** → 分辨率无关，**源/PML 伪影**，不是色散 |

### 5.2 文献

- **Oskooi et al. 2008**（PML 综述）：PML 对斜入射/倏逝波有反射，有限厚度下形成驻波；
- **Taflove Ch.7**（ABC/PML）：同理；CW 源 + 有限域 → 驻波条纹（与旧图里的斜条纹一致）；
- 本项目 CW 激光源的有限孔径/边缘衍射也会让"入射"偏离理想平面波。

**结论**：D 类是 640² 之后的主要限制；必须先修测量（RCS 减数值入射、查源/PML、近场移 3R–5R）。

---

## 6. 来源 E：Mie 级数/单位约定（本项目 `mie_cylinder.py`）

- 截断 `n_max = ceil(ka + 4ka^{1/3} + 8)`（对 ka≈6.28 足够，误差可忽略）；
- 2D 散射宽度 per-unit-length；RCS 用远场式 `(4/k)|Σ c_n e^{jnφ}|²`（脚本自检通过）；
- **total vs scattered**：`field()` 含入射、`scattered()` 不含——旧图在 PEC 内画了 `scattered()` 的
  解析延拓（物理上 E=0，**不该画**），已修正；
- 归一化：E0=1e6；By 与 Ey 的 B/E 关系（c=1 归一化）——TM 的 `|E_sc|` 由 By 的导数换算。

---

## 7. 来源 F：本项目特有（这是"为什么 C 侧还不如 c0"的直接答案）

| 机制 | 源码 | 效果 |
|---|---|---|
| `ly` 掩膜（四边任一覆盖→冻结 Ey 节点） | `MarkUpdateECellsECT`；`WarpXInitData.cpp:1679` 只传 `edge_lengths` | 切面节点被冻结 → 表面电流截断 → **m=0 低 35%**（c0 只低 14%） |
| `A_stab` 抬升（0.5·dx·dz） | `EmbeddedBoundaryInit.cpp` | 小对偶面欠松弛 → 相位偏移 + 局部幅值偏差 |
| 对偶长度一阶（scheme 2/4 实测 1.0 阶） | 同文件 | 结构残差 24%（vs c0 11%） |
| BCK 兜底（B 侧） | `WarpXFaceExtensions.cpp` | 与 c0 共有，非 C 侧特有 |
| TE 借用（G2 后） | `68a4d6f` 等 | 改变 B 侧表面响应 → TE 与 TM 的差（TE 用在 Bx/Bz，TM 用 By） |
| 无 F 项、J 项普通 | `EvolveE.cpp` | 当前算例 J=0、无 F，影响为零 |

---

## 8. 可执行清单（按收益排序）

1. **改掩膜**：`Ey` 判据从 `ly` 改为 `A_y <= A_stab`（DM 冻结；代码已备好，调用处传 `dual_A_y`）→
   恢复被冻结的切面节点，预期把 m=0 从 0.65 拉回 ≥0.85。**这是 TE 误差的第一杠杆。**
2. **解析 φ 精确 `Lx/Lz`**（+B 侧同源）：把几何从 1.0 阶提到 ≥1.8 阶（纯几何预检可先做）。
3. **测量修正**：RCS 减数值入射；真空地板的来源（源/PML）专项；近场判据移 3R–5R。
4. **相位**：在报告中把"色散相位"与"几何相位"分开（α 分解法即可）。
5. **报告口径**：所有 C 侧算例同时给"全网格 max|E|、PEC 内 max|E|、环带 L2、真空地板、α 分解"五项。

---

## 9. 参考文献

1. A. Taflove, S. Hagness, *Computational Electrodynamics: The FDTD Method*, 3rd ed., 2005（Ch.3/4/7/8）。
2. K. Yee, "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media," IEEE TAP 14(3), 1966（DOI 10.1109/TAP.1966.1138693）。
3. S. Dey, R. Mittra, "A locally conformal finite-difference time-domain (FDTD) algorithm for modeling three-dimensional perfectly conducting objects," IEEE MWCL 7(9):273–275, 1997（DOI 10.1109/7260.622821）。
4. W. Yu, R. Mittra, "A conformal finite difference time domain technique for modeling curved dielectric surfaces," IEEE MWCL 11(1):25–27, 2001（[Scholar](https://scholar.google.ro/citations?view_op=view_citation&hl=ro&user=43RLFSYAAAAJ&citation_for_view=43RLFSYAAAAJ:u5HHmVD_uO8C)）。
5. S. Benkler, N. Chavannes, N. Kuster, "A new 3-D conformal PEC FDTD scheme with user-defined geometric precision and derived stability criterion," IEEE TAP 54(6):1843–1849, 2006（DOI 10.1109/TAP.2006.875494）。
6. T. Xiao, Q. H. Liu, "Enlarged cell technique for adaptive FDTD," IEEE TAP 52/55（2004/2008）。
7. A. Farjadpour et al., "Improving accuracy by subpixel smoothing in the FDTD," Opt. Lett. 31(20):2972–2974, 2006（DOI 10.1364/OL.31.002972）。
8. A. Oskooi et al., "The failure of perfectly matched layers, and towards their redesign," J. Comput. Phys. 227:5971–5985, 2008（DOI 10.1016/j.jcp.2008.01.045）。
9. Meep 文档 `doc/docs/Subpixel_Smoothing.md` 与源码 `src/anisotropic_averaging.cpp`（本机）。
10. 本项目：`mie_cylinder.py`、`analyze_mie_convergence.py`、`EvolveE.cpp`（ECT 内核）、`EmbeddedBoundaryInit.cpp`（对偶几何）、`FieldIO.cpp`（格心平均）。

*本文是差异来源调研；实测数字来自 `/tmp/mie_now`（160² TE/TM 与真空）与 `/tmp/mie640`（640²）。*
</content>
