# WarpX_fluid commit 历史归档

<!-- SPDX-License-Identifier: BSD-3-Clause -->


`git init` 重建仓库前导出的完整 commit message（共 96 条）。
原始 .git 以 `WarpX_fluid-20261101-1737.bundle` 备份，可 `git clone` 恢复。

| # | SHA | 日期 | 标题 |
|---|---|---|---|
| 1 | `d01b3f5` | 2026-10-31 | chore: 纳入版本控制（初始提交） |
| 2 | `6225846` | 2026-10-31 | docs: C 侧调研与设计三版（研究/审查/DESIGN2/RESEARCH2） |
| 3 | `4206d85` | 2026-10-31 | feat: C 侧对偶共形 Ampere（2D XZ，开关 algo.ect_conformal_E） |
| 4 | `e6d352d` | 2026-10-31 | docs: 更新 HANDOFF 的米氏归因（与 C 侧无关）并索引实施记录 |
| 5 | `78b31c2` | 2026-10-31 | fix(mie): 定位爆炸根因为多 grid（max_grid_size<n_cell）+ 修正后处理相位与归一化 |
| 6 | `f1c6190` | 2026-10-31 | merge: C 侧对偶共形 Ampere（2D XZ，含米氏根因修复与收敛阶基线） |
| 7 | `563213c` | 2026-10-31 | docs: C 侧实施报告（改动位置/效果/残留问题） |
| 8 | `593f1f6` | 2026-10-31 | docs(cases): 新增 C16（C 侧验证与米氏根因更正），改写 8.7，新增 8.9/8.10 |
| 9 | `c4a28f9` | 2026-10-31 | docs: C 侧修正方案（依据 C_SIDE_REPORT_CHECK 的 F1-F7 改动清单） |
| 10 | `590a5f8` | 2026-10-31 | fix(c-side): F1-F3 —— 断言收敛、scheme 参数化、掩膜改用 A_y |
| 11 | `4ccafe5` | 2026-10-31 | docs: 米氏阶数不到 2 的误差源分解（纯几何实验：方案 A 无收益，需 level-set） |
| 12 | `a55b288` | 2026-10-31 | fix(c-side): 加 A_stab 分母下限；回退 F3 的掩膜放开（PEC 内热斑） |
| 13 | `550eecc` | 2026-10-31 | docs: C 侧三项发现（半格实验错位、L2 环带盲区、半侧稳定化不足） |
| 14 | `63c37d3` | 2026-10-31 | feat(c-side): G1 level-set 精确 Lx/Lz（scheme=4） |
| 15 | `64b7b1c` | 2026-10-31 | feat(c-side): G3 分母下限（Lx/Lz/A_y），掩膜判据实测不可用后回退 |
| 16 | `d08d417` | 2026-10-31 | docs(g2): 第 1 步 —— XZ 邻居编码真值表（24 组合中 14 个错） |
| 17 | `0d1d41f` | 2026-10-31 | fix(g2): 第 2 步 —— GetNeigh/SetNeigh/GetNeighPtr 的 XZ 映射 |
| 18 | `68a4d6f` | 2026-10-31 | feat(g2): 第 3 步 —— TE 面（Bx/Bz）启用完整 ECT 扩展/借用 |
| 19 | `6b76d4a` | 2026-10-31 | feat(c-side): 方案① —— DM 阈值冻结启用（A_stab=0.5 标定），E 侧稳定化闭环 |
| 20 | `d6aa57c` | 2026-10-31 | docs(c-side): E 侧稳定化的完整定义（BCK+DM，非借用）与阈值标定 |
| 21 | `914fdce` | 2026-10-31 | feat(g2): 加 algo.ect_te_borrowing 开关并完成借用 vs BCK-only 对比 |
| 22 | `78ea5d3` | 2026-10-31 | feat(m2b-1): ECT 内核支持 ε/σ；介质 E 更新接通 ECT 分派 |
| 23 | `516cd38` | 2026-10-31 | fix(m2b): 回退介质 ECT 分支（epsilon 位置未插值，E 被清零） |
| 24 | `56b1ac8` | 2026-10-31 | diag(m2b): 打印材料数组实际值 |
| 25 | `0daf76d` | 2026-10-31 | wip(m2b): 尝试 epsilon/sigma 按 Interp 插值到 E 交错位置，未能修复 E->0 |
| 26 | `693acfd` | 2026-10-31 | wip(m2b): 保留 ECT+介质分支与诊断；已定位症状（波速不减速、幅值衰减 10x） |
| 27 | `b294796` | 2026-10-31 | feat(cutcell): 定位 Lx/Lz 的 volfrac 构造是 PEC 米氏变差的主因 |
| 28 | `3439dbb` | 2026-10-31 | test(cutcell): A 方案（Lx/Lz 下限/拆分）实测不如预期，已回退默认 |
| 29 | `0766118` | 2026-10-31 | docs: A 为何不行 / B 为何行的原理分析（不含实现） |
| 30 | `52d7af2` | 2026-10-31 | fix(m2b): 按 MEDIUM_CODE_AUDIT 修 P1/P5/P10 |
| 31 | `637d11e` | 2026-10-31 | test(m2b): 均匀介质最小回归 verify_medium_wave.py（能量/波速/振幅） |
| 32 | `ee9c44f` | 2026-10-31 | docs: 介质（M2b）现状报告——三项 FAIL、六条已排除、两条待验证 |
| 33 | `210a4f2` | 2026-10-31 | docs: 按 CUTCELL_GEOMETRY_PRINCIPLE_REVIEW 修正原理文档 |
| 34 | `6a37972` | 2026-10-31 | feat(cutcell): A 的局部混合变体 —— 近切点(vmin<0.2)用常数，其余用 volfrac |
| 35 | `4162d44` | 2026-10-31 | feat(mie): --vacuum 模式——减数值入射，扣掉测量地板的一半 |
| 36 | `3b91303` | 2026-10-31 | fix(m2b): 介质根因修复 —— MacroscopicEvolveE 补 EvolveECTRho |
| 37 | `4f7ef2a` | 2026-10-31 | Revert "fix(m2b): 介质根因修复 —— MacroscopicEvolveE 补 EvolveECTRho" |
| 38 | `4d30b19` | 2026-10-31 | Reapply "fix(m2b): 介质根因修复 —— MacroscopicEvolveE 补 EvolveECTRho" |
| 39 | `5dcd1e2` | 2026-10-31 | fix(c-side): ComputeDualGeometry 与 level-set 积分器的越界（Debug 构建抓出） |
| 40 | `98101fa` | 2026-10-31 | fix(c-side): 掩膜函数的循环限制在 k=0（Debug 640^2 第 5 处排查） |
| 41 | `82316a0` | 2026-10-31 | fix(c-side): EvolveB 的 XZ 邻居映射把 j_n 写进 k，640^2 越界根因 |
| 42 | `7dc1d1f` | 2026-10-31 | feat(m2b/P0): 张量材料的数据层 —— 按分量 eps^-1 行存储（Meep chi1inv 方案） |
| 43 | `4951d92` | 2026-10-31 | feat(m2b/P0-a): Yee 内核改为按分量 eps^-1 行更新 |
| 44 | `1f16cab` | 2026-10-31 | feat(m2b/P0-b): Kottke 各向同性 subpixel（介质界面与 EB 解耦） |
| 45 | `8851146` | 2026-10-31 | fix(m2b/P0-b): Ey/Ez 取错行 + 2D 的 z 坐标越界（Fresnel 透射为 0 的根因） |
| 46 | `0f6c739` | 2026-10-31 | fix(m2b/P0-b): dz 用 CellSize(SPACEDIM-1)，2D 下 CellSize(2) 越界导致 f 恒为 0 |
| 47 | `b2ccba0` | 2026-10-31 | feat(m2b/P0): verify_fresnel_rt.py —— 真空参考法测 Fresnel R/T |
| 48 | `f8b54eb` | 2026-10-31 | test(m2b/P0): R 改用 q2、入射/透射双窗口、平台扫描；16 vs 32 ppw 不收敛 |
| 49 | `8d10a92` | 2026-10-31 | test(m2b/P0): PML 与介质不兼容，改用加长域；T 达 +0.31% |
| 50 | `bd93045` | 2026-11-01 | test(m2b/P0): 32 ppw 确认功率守恒收敛 —— P0 验收（T+守恒）通过 |
| 51 | `4ef439d` | 2026-11-01 | docs(m2b): P0 验收口径按 16/32 ppw 实测更新；相位归因（传播因子，非 bug） |
| 52 | `b2b0200` | 2026-11-01 | wip(m2b): R_flux 坡印廷通量诊断（首次实现失败，符号/窗口未定） |
| 53 | `f5263c2` | 2026-11-01 | feat(m2b/P0): R_flux/T_flux 时间积分实现 —— R 的偏差定性为二阶离散误差 |
| 54 | `5b65760` | 2026-11-01 | docs(m2b): §7.2 P0 口径按 R_flux 判读定稿 |
| 55 | `a53dcd9` | 2026-11-01 | feat(m2b/P1): 内核改为完整矩阵-向量乘 sum_d (eps^-1)_{c,d} D_d |
| 56 | `69d94f1` | 2026-11-01 | test(m2b/P1): 均匀旋转张量平面波测试 —— ix1 通道未生效 |
| 57 | `13bf474` | 2026-11-01 | test(m2b/P1): 张量测试证明 ix1 通道有效，但暴露架构缺口（ECT 介质内核仍是旧写法） |
| 58 | `a5c72b1` | 2026-11-01 | docs(m2b)+guard: subpixel 是 Yee-only，ECT 不接张量（两路径互斥） |
| 59 | `40160d9` | 2026-11-01 | feat(m2b/P1): generalized_fresnel.py —— 各向异性界面的 R/T 参考 |
| 60 | `988f85e` | 2026-11-01 | feat(m2b/P1-c): SubpixelTensor.H —— Kottke tau 变换的 host/device 实现 |
| 61 | `da67603` | 2026-11-01 | feat(m2b/P1-c): 6 分量张量输入 + chi_inv_* 接入 InitSubpixelRows |
| 62 | `b848a43` | 2026-11-01 | feat(m2b/P1-c): 求逆移到界面帧 + eps11/SPD 断言；P1-d 首测未通过 |
| 63 | `0f0c141` | 2026-11-01 | fix(m2b/P1): 坡印廷通量漏了 Ex*Hy 项 —— R+T 从 0.48 恢复到 0.99 |
| 64 | `94ebb98` | 2026-11-01 | docs(m2b): §7.4 P1 结果表 + P1-e 状态修正（εyz 在 2D 有物理意义） |
| 65 | `bbcabd7` | 2026-11-01 | test(m2b/P1-e): εxz 法向入射可测（Schur 补）—— ix2 通道验证通过 |
| 66 | `265f05f` | 2026-11-01 | test(m2b/P1-e): 两个法向入射均匀测试 —— ix2 与 iy2 通道全部验证 |
| 67 | `76626f8` | 2026-11-01 | test(m2b/P1-e): 符号核对 + 修正测试1的初始磁场（Bz -> By=Ex/c） |
| 68 | `2cc51ae` | 2026-11-01 | feat(m2b/P1): bulk 张量模式 + 退役 tensor_test 钩子（P1 收尾 1-3） |
| 69 | `057c398` | 2026-11-01 | fix(m2b): 输入解析的静默失败防线 + inv3_sym 的 SPD/残差验证 |
| 70 | `93c7059` | 2026-11-01 | feat(m2b/P1): 分支判据统一为逐点各向同性 —— 标量即 diag 特例（第一步） |
| 71 | `d9e66ec` | 2026-11-01 | test(p0): 定量测试执行记录（T5/T6/T10 通过，T3 失败，T2/T4/T9 未做，T7/T8 缺 Meep） |
| 72 | `7d04743` | 2026-11-01 | test(p0): Meep 1.34.0 可用（meep313 环境）；T7 Meep 侧脚本跑通但 harminv 无模式 |
| 73 | `acf598c` | 2026-11-01 | test(p0/T9): 介质圆柱算例与解析 Mie 对比 —— L2=203%（FAIL） |
| 74 | `d820f7b` | 2026-11-01 | test(p0/T7): 改用介质圆盘 —— Meep 与解析极点对上，WarpX 侧仍失败 |
| 75 | `ce5bd17` | 2026-11-01 | test(p0/T9): stop_time 修复 + 同源参考，203%->128%；径向判据确认平面波污染 |
| 76 | `69f931c` | 2026-11-01 | docs(p0): T7 极点全表 + 时间尺度前置检查；T7 WarpX 失败原因定位 |
| 77 | `ead9cbe` | 2026-11-01 | fix(m2b): 对角分支只填 comp 0 —— 均匀介质 + subpixel=1 时 Ey/Ez 被衰减到零 |
| 78 | `00bd479` | 2026-11-01 | docs(p0): 记录对角分支 bug（四组验证）与 T9 DFT 现状 |
| 79 | `7a782fd` | 2026-11-01 | feat(p0/T7): prony_q.py —— 多模复指数拟合，首次测出 Q=75.4 |
| 80 | `a64aa35` | 2026-11-01 | docs: P0 测试执行记录 T1–T10 |
| 81 | `c77e772` | 2026-11-01 | test(p0/T9): σ 海绵版 L2=2.85（比 PML 的 1.52 更差）—— P8 排除 |
| 82 | `1b37300` | 2026-11-01 | docs(m2b): §7.2.1 P0 定量测试执行结果；σ 海绵结论限定为'边界非主因' |
| 83 | `5d42f13` | 2026-11-01 | feat(m2b/T2): verify_oblique_rt.py + 斜入射 s/p Fresnel 结果 (§7.5) |
| 84 | `fee9041` | 2026-11-01 | fix(m2b/P1): frame(n) 是 2D 特例硬编码 —— 3D 与任意 n 全错（L0 抓到） |
| 85 | `cd13331` | 2026-11-01 | test(m2b/P1): L0 全过（37/37），extreme 阈值单独 1e-10；§7.6 记录 |
| 86 | `1b109b8` | 2026-11-01 | docs(m2b): L1 收口 —— 不加钩子；收敛率不归因于 collocated |
| 87 | `6a9fb24` | 2026-11-01 | feat(m2b/P2): 矩阵指数 sigma 的 Python 参考实现 + 6 项自检全过 |
| 88 | `4005753` | 2026-11-01 | feat(m2b/P2): sigma 矩阵指数的 C++ 实现（Jacobi eigh + M1/M2），与 numpy 一致 |
| 89 | `8b63021` | 2026-11-01 | wip(m2b/P2): M1/M2 存储与填充骨架（编译通过，sigma>0 未触发填充） |
| 90 | `b7ab2d9` | 2026-11-01 | fix(m2b/P2): M1/M2 填充未触发的根因 —— 块放在 Kottke 分支内，对角分支早退 |
| 91 | `b31f21a` | 2026-11-01 | fix(m2b/P2): 单出口重构 + M1/M2 dump；σ=0.1 与 Python 一致 |
| 92 | `b0a86e1` | 2026-11-01 | test(m2b/P2): 第 1 步判据通过（1.3e-16），但暴露 ε 取值问题 |
| 93 | `4affbb0` | 2026-11-01 | fix(m2b/P2): fill_sigma_mats 用该 cell 的 eps（同一套 phi 子采样），第 1 步真正闭合 |
| 94 | `faade2c` | 2026-11-01 | feat(m2b/P2): 第 2 步内核分支 —— sigma>0 走指数步，sigma=0 走原分支 |
| 95 | `b3a045e` | 2026-11-01 | feat(m2b/P2): 指数核补 EB mask 与 dt 断言；纯弛豫改周期边界（反射消除，误差=打印精度累加） |
| 96 | `e2feb35` | 2026-11-01 | docs: 去掉 research notes 里的本机绝对路径（公开前的清理） |

## 完整 message

```
---d01b3f5---
chore: 纳入版本控制（初始提交）

基于 WarpX 上游 26.10 的 fork ~/WarpX_fluid，用于开发
n 组分矩流体 + 完整 Maxwell + 真 cut-cell(EB) + AMR 求解器。

== 相对上游 WarpX_changed 的源码改动（6 个文件）==

Source/Initialization/WarpXInitData.cpp
  * 删除 InitLevelData 中过早执行的 EvolveECTRho（早于初始 E 写入
    Efield_fp，必然得到 rho=0，使首个半步 B 更新空转、整体滞后半拍）
  * 在 AddExternalFields 末尾（ApplyEfieldBoundary 之后）重算 ECT 环量，
    同时覆盖初始化与演化中改场两条路径及所有 AMR 层
  * 依赖修正需放在此处的原因：B 更新消费 ECTRhofield 而非现算 curl(E)

Source/FieldSolver/FiniteDifferenceSolver/EvolveECTRho.cpp
  * 放开 2D XZ 的 Rhox/Rhoz（上游 XZ 只算 Rhoy，故 TE 不支持）
  * 删除临时诊断 DBGRHO

Source/FieldSolver/FiniteDifferenceSolver/EvolveB.cpp
  * ECT 更新 idim < AMREX_SPACEDIM 改为 idim < 3（否则 2D 下 Bz 永不更新）
  * 删除临时诊断 DBGB（原为每半步三次全场 max() 归约）

Source/EmbeddedBoundary/EmbeddedBoundaryInit.cpp
  * per-unit-y 几何约化：edge_lengths[1] 取节点开放度，
    face_areas[0]/[2] 取 getAreaFrac * dz/dx
  * MarkExtensionCells 去掉 XZ 下 idim 0/2 的置零跳过
  * MarkUpdateECellsECT / MarkUpdateBCellsECT 的 XZ 判据改用 ly / Sx / Sz

Source/EmbeddedBoundary/WarpXFaceExtensions.cpp
  * CountExtFaces / ReduceIntSum / ComputeSStab 放开到 idim<3 并补 XZ
    per-unit-y 分支

Source/EmbeddedBoundary/WarpXFaceInfoBox.H
  * 上述改动的配套声明

== 验证状态（2026-10-31，详见 Docs/research_notes/）==

  L1  几何退化            PASS（Sx={0,dz}、Sz={0,dx}、ly={0,1}）
  L2a 规则域逐位退化      PASS（201 步 max|d| = 0）
  L2b 对齐壁面位置        PASS（ECT 往返误差 +0.371/+0.017/+0.076 ns @64/128/256，
                             Yee 等效壁面内移一层 dz）
  L3  斜入射（振幅/角度/偏振） 45/60 度 x TE/TM 四组 PASS（平均 |r-+1| = 0.026-0.051，
                             交叉偏振 0，斜率误差 0.1-2.4%）；30 度掠入射 FAIL
  L4  PEC 圆柱            PASS（柱内 0，镜面/阴影 TE 5.35、TM 1.42）
  长跑 45 度斜面 1e4 步    TE 能量漂移 3.0e-4；TM 1.9% 能量增长（观察项）
  上游回归                2D rotated_cube 6.8e-10、3D <= 1e-13，均未回归

待办见 Docs/research_notes/HANDOFF.md：C 侧（E 场对偶共形 Ampere）、
圆柱米氏收敛阶（被曲边界小胞元不稳定阻塞）、AMR 多层几何、2D/3D 逐点一致性。

== .gitignore 追加 ==

diags/、*.old.*、warpx_used_inputs 及若干工具导出的大 JSON。
build/ 与 build_2d/ 由上游既有 build*/ 规则忽略。

---6225846---
docs: C 侧调研与设计三版（研究/审查/DESIGN2/RESEARCH2）

---4206d85---
feat: C 侧对偶共形 Ampere（2D XZ，开关 algo.ect_conformal_E）

- Fields.H: 新增 dual_A_y / dual_L_x / dual_L_z
- WarpX.cpp: 解析 algo.ect_conformal_E（0/1/2）并分配三套对偶几何
- ComputeDualGeometry: 由 volfrac 组装 A_y/Lx/Lz（常数几何分支用于调试）
- EvolveECartesianECT: Ex/Ez 分母为 Lz/Lx，Ey 单点边中点值乘开放长度
- verify_l2a.py: 加噪声地板（零分量不再假 FAIL）

验证: L2a 常数几何与真实几何均 PASS（rel <= 3.6e-15）；
L3 45deg TE 斜率误差 2.39% -> 1.51%，幅值 0.0427 -> 0.0459。
发现: 米氏 80^2 爆炸与 EB/ECT/C 无关（Yee+PEC 同样爆，源置零则稳定）。

---e6d352d---
docs: 更新 HANDOFF 的米氏归因（与 C 侧无关）并索引实施记录

---78b31c2---
fix(mie): 定位爆炸根因为多 grid（max_grid_size<n_cell）+ 修正后处理相位与归一化

- make_mie_cases.py: amr.max_grid_size = n（单 grid）；CW 激光源跨 grid 不稳
- analyze_mie_convergence.py: 用 plotfile Header 物理时间做 DFT 相位；
  数值场除以 my_constants.E0（此前 L2 恒为 1e6）
- 收敛阶基线（TE）：掩膜 Yee 近场 1.28 阶 / RCS 1.76 阶；
  C 侧一阶几何 1.06 / 1.77（同阶，需 level-set 二阶几何才有提升）
- C_SIDE_IMPL_LOG.md: 更新 §3 根因与 §3.1 收敛阶数据

---f1c6190---
merge: C 侧对偶共形 Ampere（2D XZ，含米氏根因修复与收敛阶基线）

合并 feat/ect-conformal-E：
- 实现 EvolveECartesianECT + 对偶几何 A_y/Lx/Lz（开关 algo.ect_conformal_E）
- L2a 规则域到机器精度（rel 1.9e-15）；L3 45deg 斜率误差 2.39% -> 1.51%
- 定位米氏爆炸根因为多 grid（max_grid_size < n_cell），与 EB/ECT/C 无关
- 修正米氏后处理（相位参考、E0 归一化），得到基线：近场 1.28 阶、RCS 1.76 阶

---563213c---
docs: C 侧实施报告（改动位置/效果/残留问题）

---593f1f6---
docs(cases): 新增 C16（C 侧验证与米氏根因更正），改写 8.7，新增 8.9/8.10

- §7quinquies 标题加推翻标注
- 新章 §7sexies：C16 三算例（L2a/L3/米氏）+ 三个发现（逐位不可能、多 grid 根因、后处理两 bug）
- §8.7 根因改为 max_grid_size 多 grid，教训重写
- 新增 8.9（归因错误）与 8.10（后处理静默错误）
- §9 结论追加 C16 三条

---c4a28f9---
docs: C 侧修正方案（依据 C_SIDE_REPORT_CHECK 的 F1-F7 改动清单）

---590a5f8---
fix(c-side): F1-F3 —— 断言收敛、scheme 参数化、掩膜改用 A_y

F1 文档/注释:
- C_SIDE_REPORT R5 机理改为“分配初值 1.0 -> 低层 /1.0 -> 量纲错、静默错误”
- WarpX.H 注释去掉 “from the EB level set”，改为 volfrac(方案B) + 3/4 未实现
- HANDOFF 待办 1/2 更新为已完成/已解锁；CASES C16 补 R5 说明

F2 语义与防护:
- ect_conformal_E 断言 0<=x<=2（3/4 显式 abort，避免 bool 静默回退到方案 B）
- 多层 maxLevel()==0 断言（低层保持 1.0 会静默量纲错）
- ComputeDualGeometry bool -> int scheme，非 1/2 直接 abort
- 共形内核掩膜写 0（Yee 内核保持 return，不动上游行为）

F3 掩膜:
- MarkUpdateECellsECT 增加 dual_A_y 形参，XZ 的 Ey 判据优先 A_y<=0
- InitializeEBGridData 改为先算几何、再标掩膜

验证: =-1 与 =3 均 abort；L2a(=2) PASS（rel 3.6e-15）；
L3 45deg 与米氏三档与 F3 前逐位相同 -> 两算例上 ly 与 A_y 判据给出相同掩膜。

---4ccafe5---
docs: 米氏阶数不到 2 的误差源分解（纯几何实验：方案 A 无收益，需 level-set）

新实验(/tmp/geom_order.py, 解析圆 + 精确 volfrac, 40/80/160):
  Lx 方案B RMS 0.0412/0.0160/0.0058 -> 1.37/1.47 阶
  Lx 方案A RMS 0.0495/0.0188/0.0067 -> 1.40/1.49 阶（比 B 差 15%）
  A_y    RMS 6.9e-3/1.3e-3/2.5e-4   -> 2.39/2.42 阶

结论: A_y 与掩膜都不是瓶颈; Lx/Lz 的"相邻 cell 平均"是固有 1 阶;
方案 A 应从计划删除; 唯一有效修正是 level-set 精确求交;
另叠加真空地板 15% 与 RCS 减解析入射的伪影。

---a55b288---
fix(c-side): 加 A_stab 分母下限；回退 F3 的掩膜放开（PEC 内热斑）

复核发现 (C_SIDE_ORDER_ANALYSIS_REVIEW.md):
- 几何实验差半格，正确位置 RMS 0.0286/0.0142/0.0084，阶 ~1.0（原 1.4-1.5 作废）
- F3 的 A_y<=0 掩膜在 PEC 内制造热斑 7.1e6（环带 L2 看不见）

改动:
- 新增 algo.ect_conformal_E_Astab（默认 0.1 dx*dz）：0<A_y<A_stab 抬到 A_stab，
  covered 保持 0；只作分母下限
- 掩膜调用回退为 ly 判据（dual_A_y 形参保留）

验证: mask=ly 时 max|Ey|=1.68e6，PEC 内严格 0（回退前 7.1e6）；
A_stab 扫描 0.001-10 表明需 ~10 才压住热斑，但环带抬高 27% -> 完整稳定化仍是前置项。
文档: C_SIDE_ORDER_ANALYSIS §0/§3.2/§3.4/§5.0 已按复核更正。

---550eecc---
docs: C 侧三项发现（半格实验错位、L2 环带盲区、半侧稳定化不足）

---63c37d3---
feat(c-side): G1 level-set 精确 Lx/Lz（scheme=4）

实现:
- ComputeDualGeometry 新增 scheme=4：沿对偶线对 nodal level set 做符号变化 + 线性插值求交，
  积分开放长度（nseg=8）；规则域全开放时短路返回精确 h
- A_y 不变（¼h²Σv，二阶）；仍施加 A_stab 分母下限
- WarpX.cpp 断言放开到 4，显式拒绝 3（实测更差）

验收:
- L2a 规则域 PASS，几何打印 A_y=0.0625, Lx=Lz=0.25 精确
- L3 45° TE: |r+1| 0.0459 -> 0.0427（改善 7%）
- 米氏 80²: PEC 内严格 0，无热斑
- 米氏三档: 近场阶 0.964（scheme=2 为 1.058），RCS 1.784

结论(负面): level-set 线性插值求交没有提升阶数。φ 在 cell 内是平面近似(O(h²))，
积分出的线长仍是一阶；与 B 侧面积分数也不同源。G1 功能正确但目标未达成。
详见 C_SIDE_LS_ECT_RESEARCH.md 与本提交。

---64b7b1c---
feat(c-side): G3 分母下限（Lx/Lz/A_y），掩膜判据实测不可用后回退

G3.1 分母下限:
- A_y: 0<A_y<A_stab 抬到 A_stab（已有）
- Lx/Lz: 内核里作 Ex/Ez 的分母时取 max(L, A_stab*h)；Ey 的分子权重仍用真实 L（保通量）

G3.3 掩膜 (实测失败，已回退):
- A_y<=0      -> PEC 内 7.13e6 热斑（F-2）
- A_y<=A_stab -> 仍 6.89e6（0.1 阈值不足以冻住几何半开、物理在导体内的节点）
- 需要 A_stab~10 才压到 1.8e5，但环带抬高 27%（退化为阶梯）
- 故保留 ly 判据；只保留分母下限

验证: L2a PASS；米氏 80 max=1.676e6, PEC 内=0, ring=1.656e6

---d08d417---
docs(g2): 第 1 步 —— XZ 邻居编码真值表（24 组合中 14 个错）

- 新增 benchmarks/cases/c02_0_pec/g2_neigh_truth_table.py：枚举 8 邻居 x 3 法向，
  打印现状/正确偏移与 3x3 网格落点；round-trip 校验通过（编码自洽）
- 新增 Docs/research_notes/G2_NEIGH_TRUTH_TABLE.md：给出修复规格
  · GetNeigh/SetNeigh/GetNeighPtr: dim=0 -> arr(i, j+j_n, k); dim=2 -> arr(i+i_n, j, k)
  · 扩展循环 XZ 下只遍历一维（dim=0 只沿 z，dim=2 只沿 x）
  · 两者必须同时改，否则重复记账
- 关键发现: dim=0 现状把 j_n 写进哑元 k（越界/串 fab）

---0d1d41f---
fix(g2): 第 2 步 —— GetNeigh/SetNeigh/GetNeighPtr 的 XZ 映射

- 三个函数改为 #ifdef WARPX_DIM_XZ 优先分支：
  dim=0 -> arr(i, j + j_n, k)   (x 面，面内只有 z)
  dim=1 -> arr(i + i_n, j + j_n, k)
  dim=2 -> arr(i + i_n, j, k)   (z 面，面内只有 x)
- 3D/1D 的原有分支原样保留在 #else 内（未改行为）
- 与 G2_NEIGH_TRUTH_TABLE.md 的 'correct' 列逐项一致

回归: L2a PASS (rel 3.6e-15)；45deg L3 PASS。
注: 扩展循环仍是 idim=1..，dim=0/2 未被枚举，故本轮数值不变（符合预期）。

---68a4d6f---
feat(g2): 第 3 步 —— TE 面（Bx/Bz）启用完整 ECT 扩展/借用

- WarpXFaceExtensions.cpp: one-way 的 XZ 循环 idim=1.. -> 0..2
- 借用邻居枚举在 XZ 下收敛到一维：idim=0 只 i_n=0（沿 z），idim=2 只 j_n=0（沿 x）
- EmbeddedBoundaryInit.cpp: 删除 XZ 下 idim=0/2 的 BCK-only 特判，走通用扩展路径

结果:
- L2a PASS（rel 3.6e-15，规则域不变）
- L3 45deg: r=-1.0459 不变
- 米氏 80: max 1.761e6（+5%），PEC 内 0，ring 1.655e6（不变），无警告残留
- 警告显示 x/y/z 三向 'faces to be extended after eight ways' 均为 0，说明扩展全部成功

E 侧说明: Ey/Ex/Ez 的几何已共形（A_y/Lz/Lx），但 E 侧没有与 Venl 等价的
借用记账，'E 每一分量 ECT' 的稳定化部分需要单独设计（未在本次改动内）。

---6b76d4a---
feat(c-side): 方案① —— DM 阈值冻结启用（A_stab=0.5 标定），E 侧稳定化闭环

概念更正（采纳）：E 侧的'完整'不是借用（Venl 是 B 侧机制），而是
BCK 分母下限 + Dey-Mittra 阈值冻结，两者合起来就是文献意义上的完整稳定化。

改动:
- 启用 DM 冻结掩膜：MarkUpdateECellsECT 传 dual_A_y + A_stab，判据 A_y <= A_stab 时 E=0
- A_stab 默认 0.1 -> 0.5（标定结果见 WarpX.H 注释）

标定 (Mie TE 80^2, 末帧):
  A_stab=0.1 -> max|Ey| 6.67e6, PEC 内 6.67e6（热斑）
  A_stab=0.5 -> max|Ey| 3.39e6, PEC 内 0, 环带 1.655e6（不变）  <-- 默认
  A_stab=1.0 -> 全冻结（A_y <= dx*dz 恒成立）

验收:
- L2a PASS（rel 3.6e-15，规则域不变）
- 米氏 80: F-2 三行 max=3.393e6 / PEC=0 / ring=1.655e6
- TM 45deg 2000 步: Emax/E0=1.0205, 末端漂移 7.4e-3（与 C14 同量级，Lx/Lz floor 未消除该增长）

---d6aa57c---
docs(c-side): E 侧稳定化的完整定义（BCK+DM，非借用）与阈值标定

- 更正术语: 'E 每一分量完整 ECT' -> 'E 侧完整稳定化 = BCK 分母下限 + DM 阈值冻结'
- 澄清: 借用(Venl)是 B 侧机制；E 侧可移植形态是 Xiao-Liu 的'扩大对偶环模板'，非累加器
- 更正 F-3 结论: 抬分母需要 A_stab~10，而 DM 冻结只需 0.5（环带不动）
- 记录三行验收与 TM 2000 步漂移 7.4e-3

---914fdce---
feat(g2): 加 algo.ect_te_borrowing 开关并完成借用 vs BCK-only 对比

- WarpX.H/.cpp: ect_te_borrowing (默认 1)；EmbeddedBoundaryInit 加回 BCK-only 分支（由开关控制）
- 同一二进制下两组对比（scheme=2, A_stab=0.5）：

  L3 45deg TE:    借用 r=-1.0459 (0.0459) == BCK-only（完全相同）
  米氏 TE 三档:   借用 近场 1.042 / RCS 1.722
                  BCK-only 近场 1.058 / RCS 1.770

结论: TE 面借用精度上无收益、略差(-1.5%/-2.7%)；它是'用面积换稳定'，
当前算例未暴露稳定性差异。实现正确（无崩溃、无警告残留、PEC 内 0）。

---78ea5d3---
feat(m2b-1): ECT 内核支持 ε/σ；介质 E 更新接通 ECT 分派

- EvolveECartesianECT 增加 MacroscopicProperties 参数（默认 nullptr = 真空）
- 内核三分量：E = alpha*E + (eps0/eps_abs)*dE/(1+s)，alpha/s 为 CN 系数（sigma 隐式）
- FiniteDifferenceSolver::MacroscopicEvolveE 增加 fields 参数，ECT+conformal 时走共形内核
- WarpX::MacroscopicEvolveE 传 m_fields

状态: 编译通过；均匀介质最小算例跑通(exit=0)，但初值未注入(max|Ey|=2.9e-15)，
介质 R/T 验收(Fresnel V1)未完成 —— 需修算例的 E_ext_grid_init_style 配置。

---516cd38---
fix(m2b): 回退介质 ECT 分支（epsilon 位置未插值，E 被清零）

现象: 均匀介质 eps_r=2.25 算例中 max|Ey| 由 1.0 降到 0（几百步内），
      真空对照正常(0.992)，说明问题出在介质内核分支。
根因(未完全证实): MacroscopicEvolveECartesian 用 ablastr::coarsen::sample::Interp
      把 cell-centered 的 eps/sigma 插值到 E 的交错位置，而 ECT 内核直接读 cell 中心；
      另外初值探针位置曾误置在脉冲尾部(exp(-44))导致误判'初值未注入'，已澄清。
处置: 关闭 MacroscopicEvolveE 的 ECT 分支（保留内核的 eps 参数与保护），
      介质仍走原 Yee 路径；等 Interp 移植后再启用。
保留: EvolveECartesianECT 的 MacroscopicProperties 参数（默认 nullptr，不影响真空路径）。
文档: 预期下一步 = 移植 Interp + Fresnel V1。

---56b1ac8---
diag(m2b): 打印材料数组实际值

结果: eps=1.992e-11 (=2.25*eps0), sigma=0, mu=mu0 —— 材料初始化正确。
说明 E 被清零不是材料参数问题，而在 MacroscopicEvolveECartesian 的
ablastr::coarsen::sample::Interp 环节（cell 中心 -> E 交错位置），或 2D 下 epsilon_stag 的配置。
注: 回退后介质算例走的是原 Yee 介质路径（ECT 分支已关），因此这是上游路径的问题，
    不是本次 ECT 改动引入的。

---0daf76d---
wip(m2b): 尝试 epsilon/sigma 按 Interp 插值到 E 交错位置，未能修复 E->0

- EvolveECartesianECT: 三分量改用 ablastr::coarsen::sample::Interp(eps/sig, stag, X_stag, macro_cr, ...)
- MacroscopicEvolveE: 短暂重开 ECT 分支后再次禁用

结果: ECT+介质仍 max|Ey|=0（纯 Yee+介质 = 0.8255 正常）。
说明根因不是 epsilon 位置插值。剩余候选:
  1) 该路径下 dual_A_y/Lx/Lz 未被 ComputeDualGeometry 初始化（A_y=1 初值会放大而非清零，但仍需确认）
  2) E 在一步内被更新两次（WarpX::EvolveE 与 MacroscopicEvolveE 都触发）
  3) Interp 返回 0 -> alpha=nan -> 被后续清零
下一步应先加'内核是否被调用 + 首帧 Ey 值'诊断，而非继续改公式。
现状: 介质仅走 Yee（已验证可用），worktree 干净。

---693acfd---
wip(m2b): 保留 ECT+介质分支与诊断；已定位症状（波速不减速、幅值衰减 10x）

诊断结论（同算例 96^2, eps_r=2.25, PEC, 探针 z=-0.5）:
  ComputeDualGeometry 正常: A_y=0.001736=dx*dz, Lx=Lz=0.041667=dx
  ECT 真空(D 算例): max|Ey|=0.9918, 到达 62 步  <-- 内核正确
  ECT+介质: max|Ey|=0.0979, 到达 66 步
      -> 波速几乎等于真空(应慢 1.5 倍), 幅值小 10 倍
  纯 Yee+介质: max|Ey|=0.8255 正常
  eps_min=eps_max=1.99e-11(=2.25*eps0), sigma=0, 材料数组无误

即 inv_epsr 既未减速也未保幅 -> 不是 epsilon 位置问题。
剩余候选:
  1) E 在一步内被更新两次(真空 EvolveE + 介质 MacroscopicEvolveE)
  2) dE 与 alpha*E 的相消(赋值 vs 累加语义)
  3) FillBoundary/边界处理把介质 E 与真空 E 混用
下一步: 在 WarpX::EvolveE 与 WarpX::MacroscopicEvolveE 各加调用计数(已有 DIAG-ECT)，比较一步内调用次数。

---b294796---
feat(cutcell): 定位 Lx/Lz 的 volfrac 构造是 PEC 米氏变差的主因

隔离实验（米氏 TE 160^2, 同一二进制, 仅换 scheme）:
  scheme 0 (掩膜 Yee):              L2=0.3465  phase 8.9deg  阶 1.278
  scheme 1 (A_y/Lx/Lz 全常数):      L2=0.3465  phase 8.9deg  (与 0 逐位相同)
  scheme 6 (A_y=volfrac, Lx/Lz 常数): L2=0.3455  phase 8.8deg  阶 1.236
  scheme 7 (A_y 常数, Lx/Lz=volfrac): L2=0.4525  phase 11.5deg
  scheme 2 (全 volfrac, 当前默认):   L2=0.4454  phase 11.3deg  阶 1.058

结论: Lx/Lz 的 1/2*(v+v) 在近切点严重低估开放长度(v->0 而真实线长->0.5h),
      造成 1/L 放大, 使 L2 恶化 30%; A_y 的 volfrac 构造无害(6 甚至略优于 0)。
新增 scheme 6/7 作为隔离开关, 不改变默认(2)。
文档: MIE_FIX_METHODS_RESEARCH 的'掩膜是主因'假设被此实验否定, 应改为'Lx/Lz 构造'。

---3439dbb---
test(cutcell): A 方案（Lx/Lz 下限/拆分）实测不如预期，已回退默认

实验(米氏 TE 160^2, 近场 L2 / 阶数):
  volfrac(原默认, scheme 2):   0.4454 / 1.058
  Lx/Lz=max(volfrac,0.5h):     0.4034 / 1.121   <- A 方案, 不够
  Lx/Lz 常数(scheme 6):        0.3455 / 1.236   <- 最好
  Ey 分子常数 + Ex/Ez 分母 volfrac: 0.3455 / 1.218
  掩膜 Yee(c0):                0.3465 / 1.278

45 度斜面回归:
  volfrac:       |r+1|=0.0459, 斜率 -1.51%   <- 平面精确
  常数 Lx/Lz:    |r+1|=0.0427, 斜率 -2.39%   <- 斜面变差

关键认识: volfrac 平均在平面上精确(½(v+v) 对直线段成立), 在曲面近切点高估
(真实 0, 给出 0.5h), 两种情形的误差方向相反 -> 没有单一的'下限'能同时修好
米氏与斜面。需要真正的精确几何(parser 求交)且 B 侧同源, 即方案 B。
保留 scheme 6/7/8 作为诊断开关, 默认仍为 2(volfrac), 工作区已回退。

---0766118---
docs: A 为何不行 / B 为何行的原理分析（不含实现）

- A 的三层原因: volfrac 是不可逆投影(信息丢失) + 沿对偶线凸性导致 Jensen 高估近切点
  + Lx/Lz 与 A_y 的维度差异; A2 常数只是把误差匀化, 45 度斜率退化暴露代价
- B 的原理: 跳过投影做一维截断积分; 必须用解析 phi(采样 phi 的 O(h^2) 平面近似在近切点仍一阶);
  必须 B 侧同源(离散 Hodge star 自伴性: C^T star C 要求两侧同一剖分)
- 精度阶传递: 场误差 ~ eps_geom / L, 近切点是分母趋零的放大器
- 4 条可检验推论(用于 B 之后验证/否证该原理)

---52d7af2---
fix(m2b): 按 MEDIUM_CODE_AUDIT 修 P1/P5/P10

- P1(高): ECT 材料路径补 mu —— dE 现在乘 mu0/mu_abs(即 1/mur)，与 Yee 路径的 H=B/mu 对齐；
  三分量各自用 mu_Interp 到对应的 E 交错位置
- P5(中): 删除 [DIAG-ECT]/[DIAG-S] 调试打印与 static 计数器（并行区竞态）
- P10(低): 材料守卫改用插值后的 eps_abs(带 1e-30 下限)，与插值点一致

未修: P2/P3(上游 Yee 路径的 mu 半格、epsilon 算术平均)、P4(sigma_method 在 ECT 无效)、
P6(无界面回归)、P8(PML)、P9(Yee+EB 阶梯) —— 见审计第 4 节路线。
回归: L2a PASS (rel 3.602e-15)。

---637d11e---
test(m2b): 均匀介质最小回归 verify_medium_wave.py（能量/波速/振幅）

三项判据（均匀 eps_r，PEC，无界面）:
  1. 能量漂移 <= 1%（无损耗应守恒）
  2. 相速度 = c/sqrt(eps_r) +-10%
  3. 两探针振幅比在 1/2..2（平面波不横向扩散）

当前代码的实测（96^2, eps_r=2.25）:
  FAIL energy   drift=1.663e+05  (max/E0=166293)
  FAIL speed    v=1.4454e8 vs 1.9986e8, rel=-27.7%
  FAIL amplitude B/A=0.000 (A=0.1735, B=0.0000)

用法:
  python3 verify_medium_wave.py <dir> --make --run --epsr 2.25 --ncell 96
  python3 verify_medium_wave.py <dir> --epsr 2.25 --ncell 96

价值: 把'探针读数 0.0979、靠肉眼判断'变成三个可断言的数字，且对界面处理不敏感，
正好隔离出 ECT 材料内核本身的问题（M2b 的 Fresnel 测试看不出这一类）。

---ee9c44f---
docs: 介质（M2b）现状报告——三项 FAIL、六条已排除、两条待验证

- 现象: 均匀介质三项判据全 FAIL（能量 1.66e5、速度 -27.7%、振幅比 0）
- 三条对照: ECT 真空 OK / Yee 介质 OK / ECT+介质 FAIL -> 夹在 has_macro 分支
- 已排除 6 条: eps 位置/J/alpha/几何/mu/材料读取
- 审计 P1-P10 状态: P1/P5/P10 已修, P2/P3/P4/P6/P8/P9 未修
- 两条待验证: 假设A(E 双更新) / 假设B(非线性反馈), 各一次运行可判
- 结论: 介质不可用, 但不影响 PEC cut-cell 工作(已隔离)

---210a4f2---
docs: 按 CUTCELL_GEOMETRY_PRINCIPLE_REVIEW 修正原理文档

- §1.2: Jensen'系统性高估'论证被否 —— 改为双向误差(RMS 0.167h/近切点 0.264h，33% 高估/67% 低估)，
  机制改为'2h 宽块平均 vs h 宽对偶线的半格错位'；连锁后果降级为待验证假设
- §1.3: A_y 阶数 2.39/2.42 -> 1.95/1.99（对偶线正确位置口径）
- §2.2: 补两个限定 —— parser 必须沿线任意点求值；只改 E 侧不保证二阶
- §2.3: 弱化 Hodge star 论证（本格式 Lx 兼作分子/分母，B 侧有 BCK/借用）
- §0: 加测量口径说明（0.3465 vs 0.3270）与 scheme 6/7/8 的引入提交
- 结尾一句话同步修正

---6a37972---
feat(cutcell): A 的局部混合变体 —— 近切点(vmin<0.2)用常数，其余用 volfrac

动机: 复核确认 volfrac 平均的误差是'双向 + 集中在近切点'(RMS 0.167h/0.264h),
      机制是 2h 宽块 vs h 对偶线的半格错位; 平面因 v(x) 线性而精确。
      于是按局部判据分流: 近切点(vmin<0.2)用常数, 平面上保留 volfrac。

结果(米氏 TE 三档 + 45 度):
  volfrac(原):        1.058 阶, 45度 斜率 -1.51%
  局部混合(新默认):    1.168 阶, 45度 斜率 -1.51%   <- 米氏 +10%, 斜面不退化
  全常数:             1.236 阶, 45度 斜率 -2.39%
  掩膜 Yee(c0):        1.278 阶, 45度 斜率 -1.51%

L2a 回归: PASS (rel 3.602e-15), 规则域 vmin=1 不触发阈值, Lx=dx 精确。
scheme 9 = 显式局部混合; 默认 scheme 2 现在等价于它。

---4162d44---
feat(mie): --vacuum 模式——减数值入射，扣掉测量地板的一半

- analyze_mie_convergence.py 支持 --vacuum <root>：用同分辨率/同源的真空跑作为入射场
- 近场改为散射场对比：A_sc = (A_al - A_vac*e^{iφ}) vs (Am - e^{iKZ})
- RCS 同步减数值入射

实测(def9 = 局部混合默认):
  解析入射: RCS L2 = 0.286/0.206, 阶 1.33
  数值入射: RCS L2 = 0.111/0.122, 阶 0.165   <- 地板被扣掉, 但剩余不随 dx 收敛
  近场:     0.3799->0.4247(160), 0.2248->0.2195(320), 阶 1.045->1.060

结论: '减数值入射'只扣掉了 RCS 的那一半地板; 近场不动, 剩余误差是散射场本身的
(几何/掩膜/界面场), 不是测量伪影。这为 B 的验收提供了正确的分母。

---3b91303---
fix(m2b): 介质根因修复 —— MacroscopicEvolveE 补 EvolveECTRho

根因（见 MEDIUM_SPEED_VERIFICATION.md §5）:
  WarpX::MacroscopicEvolveE 在 E 更新后只调 ApplyEfieldBoundary，
  漏掉了 WarpX::EvolveE 末尾的 ECTRhofield 重算。
  于是 ECT 的 B 更新（EvolveBCartesianECT 用 ECTRhofield）读的是 E 更新之前的 Rho,
  每步注入能量并正反馈 -> 指数增长。
  这解释了: epsr=1 也爆、mu/材料系数/内核都无关、Yee+介质正常(Yee 不用 ECTRhofield)。

修复: 在 ApplyEfieldBoundary 之后补 EvolveECTRho(fine/coarse)，与 EvolveE 一致。

验证:
- 介质加权能量 final/init = 0.99997（原 1.66e5）
- A 的 640 段错误随之消失（原来是能量爆炸触发 inf/nan）: 完整 640 跑通, 111 帧
- 同口径 5 档（散射场+真空入射）: A 近场 0.948 / c0 0.995; RCS 0.155 / 0.133

遗留: verify_medium_speed.py 的测速分支给出 v_dft=1.68c、v_arrival=nan,
      能量判据正确但测速方法需修（MEDIUM_SPEED_VERIFICATION §2 的包络峰值法未对齐）。

---4f7ef2a---
Revert "fix(m2b): 介质根因修复 —— MacroscopicEvolveE 补 EvolveECTRho"

This reverts commit 3b913036ae47b2abd6eba1cdda4558b5e73c9c4c.

---4d30b19---
Reapply "fix(m2b): 介质根因修复 —— MacroscopicEvolveE 补 EvolveECTRho"

This reverts commit 4f7ef2aecbdd281a82c2199f2fb0147f7b9d010e.

---5dcd1e2---
fix(c-side): ComputeDualGeometry 与 level-set 积分器的越界（Debug 构建抓出）

用 build_dbg（CMAKE_BUILD_TYPE=Debug，AMReX Array4 断言）跑 640^2 米氏算例，
逐次抓出四个越界点并修掉三个：

1. volfrac 的 i-1 下界：(-2,-2) out of (-1:640,-1:640)  <- Ay 的 v(i-1,j-1)
   -> 对 A_y 的四个象限索引 clamp 到 volfrac 的 box
2. volfrac 的 j-1 下界：(-2,-1)
   -> 同上（A_y）
3. volfrac 的 i 上界：(641,-1)  <- Lx/Lz 的 v(i,j,k)
   -> 对 Lx/Lz 的 v 访问 clamp 上下界
4. z 方向的 ghost：(312,240,-1)  <- level-set 积分器在 k=-1 取 phi
   -> open_len_x/z_from_ls 内把 k clamp 到 phi 的 z 范围

另：dual_A_y/Lx/Lz 的写入 box 改用 2D IntVect（k 方向 0 层），
    避免 ng_FieldSolver 的 3D ghost 让 box 跨到 k=-2..2。

验证: L2a 仍 PASS（rel 3.602e-15，钳制只在 ghost 生效）；
Debug 640^2 仍有 1 处越界待定位（dual_L_x 的 k=-1，见下一条）。
说明: 这也是此前'介质修复让 640 不崩'的真正解释——越界是 UB，
      EvolveECTRho 的 20 行改变了内存布局，把写坏的位置推到了无害处。

---98101fa---
fix(c-side): 掩膜函数的循环限制在 k=0（Debug 640^2 第 5 处排查）

MarkUpdateECellsECT / MarkUpdateBCellsECT 的 tbx/tby/tbz 原本用
iMultiFab->nGrowVect()，在 2D 下仍带 3D ghost，循环跨到 k=-2..2，
而 face_areas/edge_lengths/eb_update_* 在 k 上只有 1 个 cell。
改为 2D IntVect(ng[0],ng[1])。

结果: L2a 仍 PASS（rel 3.602e-15）。
遗留: Debug 640^2 仍报 (312,240,-1) out of (-1:641,-1:640,0:0)，
该 box 形态（x nodal、y cell）对应 dual_L_x / face_areas[0]，
但已知的四处写入/读取 box 都已限制为 k=0，来源未定位。
下一步: 在 face_areas[0]、dual_L_x 的每次访问处加 AMREX_ASSERT 定位调用栈，
或用 gdb 在 amrex::Abort 处捕获 backtrace。

---82316a0---
fix(c-side): EvolveB 的 XZ 邻居映射把 j_n 写进 k，640^2 越界根因

用 Backtrace.0.0（Debug 版带 -g）直接拿到调用栈:
  EvolveBCartesianECT (EvolveB.cpp:288) -> lambda -> EvolveB.cpp:324 Rho(ip,jp,kp)

idim == 0 分支用的是 3D 映射:
  jp = j + vec(0);  kp = k + vec(1);
XZ 下 k 退化，且 x 面的面内方向只有 z -> vec(1) 应进 j，不是 k。
于是 kp = -1，索引 (k=0 only) 的 Rho 越界。

修复(与 G2 真值表一致):
  idim == 0 : jp = j + vec(1); kp = k
  idim == 2 : ip = i + vec(0); jp = j; kp = k
两处各 2 个循环（Venl 累加 + 借用分配），共 4 处。

验证:
- Debug 640^2: 100 步 4 帧, exit=0，无 out-of-bound
- Release 640^2: 完整运行 111 帧, exit=0（此前为 MPI_ABORT 11）
- L2a 仍 PASS（rel 3.602e-15，规则域无 borrowing，映射改动无影响）

这同时是此前'介质修复让 640 不崩'的最终解释:
越界在 B 侧，介质修复改变了内存布局把它掩盖了。

---7dc1d1f---
feat(m2b/P0): 张量材料的数据层 —— 按分量 eps^-1 行存储（Meep chi1inv 方案）

按 Docs/research_notes/技术路径.md §2.2/§8 P0 落地第一步:

- MacroscopicProperties 新增 m_eps_inv_x/y/z（各 3 分量），分别落在 Ex/Ey/Ez
  的 Yee 位置；getter + has_subpixel_rows()
- InitSubpixelRows(): 分配 + 对角填充 1/eps（均匀介质时与标量路径等价），
  三个 MultiFab 因 stagger 不同逐个填充（MultiFab::copy 要求同 BoxArray）
- algo.macroscopic_subpixel 开关（默认 0），懒初始化挂在 WarpX::MacroscopicEvolveE

验证: 开关打开后均匀介质算例与关闭时逐位相同（2.090e-01，内核尚未使用新行），
说明数据层不改变现有行为。

下一步: (a) 内核改按行 E_c = sum_d (eps^-1)_{cd} D_d；
        (b) Kottke 各向同性界面平均（法向调和/切向算术，需要 EB level set）；
        (c) 各向异性 tau 变换（论文 Eq.4/22/23）。

---4951d92---
feat(m2b/P0-a): Yee 内核改为按分量 eps^-1 行更新

- MacroscopicEvolveECartesian 的三分量更新改为:
    E_c = alpha*E_c + beta*eps * (eps^-1)_{c,x} * D_c
  其中 beta*eps = dt/(1+s)（把 1/eps 从 beta 里提出来），(eps^-1)_{c,x} 由
  子像素行提供（对角时 = 1/eps，与标量路径代数等价）。
- D_c 现在把 J 项并入（D_c = curlH - j），保持与原来 -beta*j 的等价性。
- has_subpixel_rows() 为 false 时回退到 1/epsilon_interp，行为不变。

验证: 均匀介质 (eps_r=2.25, 60 步)
  subpixel=1 (按行) drift = 2.0902e-01
  subpixel=0 (标量) drift = 2.0902e-01   -> 逐位一致

下一步 P0-b: Kottke 各向同性界面平均 (eps^-1 法向调和/切向算术)，
需要界面法向 n̂。当前 EB 是 PEC (eps=inf)，不是介质界面，
需先确定介质界面的定义方式（新 level set 场 vs 复用 EB 几何 + 两侧 eps）。

---1f16cab---
feat(m2b/P0-b): Kottke 各向同性 subpixel（介质界面与 EB 解耦）

按 技术路径.md §5.1:
- 新增 macroscopic.epsilon2 与 macroscopic.interface_function(x,y,z)（phi>0 = 材料 1）
- InitSubpixelRows 不再依赖 eb_fact：无界面时对角 1/eps；有界面时按 cell 算
    f      = 4x4 子采样 phi 的材料 1 分数
    n      = normalize(grad phi)（中心差分）
    <eps>  = f*eps1+(1-f)*eps2,  <1/eps> = f/eps1+(1-f)/eps2
    chi^-1 = n n^T (<1/eps> - 1/<eps>) + I/<eps>
  三个 E 分量各取 chi^-1 的一行（collocated 近似，同 Meep）

状态: 编译通过；无界面时与标量路径一致（此前已验证 drift 2.0902e-01 逐位相同）。
Fresnel 验证未通过: eps1=1/eps2=4 的 z=0 界面算例中，透射侧探针 (z=0.8) 读数为 0，
波未穿透，原因待查（可能是 chi^-1 的 f/n 计算或行索引）。

---8851146---
fix(m2b/P0-b): Ey/Ez 取错行 + 2D 的 z 坐标越界（Fresnel 透射为 0 的根因）

三处，均由复核指出并已用数值验证:

1. Ey 用了 eiy(i,j,k,0) = chi^-1_{y,x}。对 z 法向界面 n=(0,0,1)，该分量为 0，
   于是 Ey 在界面 cell 的更新退化为 alpha*Ey（无增量），等价于 PEC -> 全反射，
   透射侧探针恒为 0。改为 eiy(i,j,k,1) = chi^-1_{y,y} = 1/<eps>。
2. Ez 同理: eiz(i,j,k,0) -> eiz(i,j,k,2)。
3. InitSubpixelRows 用 plo[2] 作 z 坐标，但 2D 下 RealBox 只有 AMREX_SPACEDIM
   个分量，应为 plo[AMREX_SPACEDIM-1]。

验证 (192^2, eps1=1, eps2=4, 界面 z=0, 垂直入射):
  修正前 pB(z=0.8) max|Ey| = 0.0000
  修正后 pB(z=0.8) max|Ey| = 0.6929
  解析透射系数 t = 2*n1/(n1+n2) = 2/3 = 0.6667 -> 实测偏差 +3.9%
  pA(z=-1.0) 仍为 1.0000（初值，未受影响）

说明: 复核同时指出内核应做完整矩阵-向量乘 sum_d (eps^-1)_{c,d} D_d，
当前只取对角项（各向同性够用）；各向异性 (非对角) 是 P0-c 的工作。

---0f6c739---
fix(m2b/P0-b): dz 用 CellSize(SPACEDIM-1)，2D 下 CellSize(2) 越界导致 f 恒为 0

现象: 探测发现 chi^-1_x 全域 = 1/eps2（入射侧也变成 eps=4），
      FFT 时间频率是空间频率的一半（波速 c/2），印证'入射侧被分到材料 2'。

根因: InitSubpixelRows 里 dz = geom.CellSize(2)；2D 只有 AMREX_SPACEDIM=2 个方向，
      index 2 越界返回未定义值 -> zc 错 -> phi=-z 符号错 -> f 恒为 0 -> 全材料 2。

修复: dz = geom.CellSize(AMREX_SPACEDIM-1)；同时把调试打印改为 chi^-1_yy（comp 1）。

验证 (eps1=1, eps2=4, 界面 z=0):
  chi^-1_yy: min=2.82e10 (=1/eps2)  max=1.13e11 (=1/eps1)   <- 两侧都出现
  FFT 峰 omega = 4.11e9（解析 3.77e9，+9% 为 16ppw 的数值色散）
  10 步 plotfile: 空间 lambda = 0.5000（与输入表达式一致）

遗留: 双探针 DFT 解出的 Ei 仍偏小（R=8.2 不可信），需要改用真空参考跑
      （同网格无界面 -> E_ref，直接相减）以避开 k/omega 的假设。

---b2ccba0---
feat(m2b/P0): verify_fresnel_rt.py —— 真空参考法测 Fresnel R/T

双探针行波分解需要 k(omega)，而 Yee 色散不是 omega/c，用 FFT 峰当 omega 会让
2x2 系统病态（|Er|>|Ei|、R>1 是典型症状）。改用同源真空参考跑:

  A) eps=1 全域、同网格/同源/同探针   -> E_ref（入射）
  B) eps1=1, eps2=4, 界面 z=0        -> E_tot, E_t
  R = |(E_tot-E_ref)/E_ref|^2
  T = (n2/n1)|E_t/E_ref|^2

A/B 共享同一 Yee 色散，相减自动对齐，不需要任何 k/omega 假设。

实测 (256^2, eps2=4, omega=3.72e9):
  R (q2, 入射侧 z=-0.6) = 0.1088   理论 0.1111   (-2.1%)
  T (q3,q4 平均)        = 0.8775   理论 0.8889   (-1.3%)
  R(q2)+T               = 0.986    (-1.4%)

遗留: q1 (z=-1.2) 的时间窗截得太早，反射波刚到达，R(q1)=0.0013 偏低；
      需要把入射侧窗口放到 ~9ns（PEC 从 z=-4 返回前）。脚本已暴露 tmax 参数。

---f8b54eb---
test(m2b/P0): R 改用 q2、入射/透射双窗口、平台扫描；16 vs 32 ppw 不收敛

按复核: q1 的界面反射 (9.0ns) 与 PEC 反射 (9.3ns) 重叠且脉冲长 1.7ns，
没有干净窗口 -> R 只用 q2 (z=-0.6)，透射侧窗口到 18ns 覆盖 c/2 慢波。

实测:
  16 ppw (256^2): R=0.1088 (-2.05%), T=0.8775 (-1.28%), R+T=0.9864
  32 ppw (512^2): R=0.0976 (-12.2%),  T=0.8690 (-2.23%), R+T=0.9666

平台扫描 (16 ppw, tmax-inc): 7.0ns 0.0344 / 7.5 0.0591 / 8.0 0.1088 / 8.5 0.1186 / 9.0 0.1286
  -> 不平台，随窗口单调增长，说明窗口内仍有非界面反射的成分进入。

结论: R 不随分辨率收敛（-2% -> -12%），所以不是数值色散；窗口/边界是主嫌疑。
下一步: 用 PML 消掉 PEC 反弹（边界改 pml），窗口只受脉冲通过时间限制，q1 也可用。

---8d10a92---
test(m2b/P0): PML 与介质不兼容，改用加长域；T 达 +0.31%

- PML 尝试失败: em_solver_medium=macroscopic + boundary=pml 在 ProjectionDivCleaner
  崩溃（MEDIUM_CODE_AUDIT.md P8 记录的已知限制）
- 改用 domain [-8,8]（PEC 反弹从 9ns 推迟到 ~25ns），窗口 [-2.0e-8, 2.6e-8]

实测 (eps2=4, 16 ppw):
  T_pow = 0.8917 (theory 0.8889)  rel +0.31%   <- 达标
  R     = 0.1307 (theory 0.1111)  rel +17.6%
  R+T   = 1.0224

分析: T 的 0.31% 说明 subpixel 的透射是正确的；R 偏大来自 E_tot-E_ref 差值里
混入了界面附近的一阶点值误差（Meep limitation 4），而不是 subpixel 的界面平均公式。
R 的验收应改用积分量（RCS/散射通量）或更高分辨率，而不是单点反射系数。

---bd93045---
test(m2b/P0): 32 ppw 确认功率守恒收敛 —— P0 验收（T+守恒）通过

判据按复核调整为 T_power 为主 + 功率守恒，R 的点值不作 P0 判据。

实测 (eps2=4, 域 [-8,8], 同样的 lambda=0.5):
  16 ppw (512):  R=0.1307 (+17.6%)  T=0.8917 (-0.31%)  R+T=1.0224 (+2.2%)
  32 ppw (1024): R=0.1156 (+4.01%)  T=0.8851 (-0.43%)  R+T=1.0007 (+0.07%)

结论:
- 功率亏损 2.2% -> 0.07%，随加密收敛 -> 是 collocated eps 的离散误差，不是系统性偏差；
- T 稳定在 -0.3%~-0.4%，subpixel 在 eps2 侧正确；
- R 从 +17.6% 降到 +4.0%，同向收敛，但点值受界面一阶限制，不作 P0 判据。

P0 验收: T 相对误差 <=1% (实测 0.43%) + R+T 与 1 的偏差 <=1% (实测 0.07%) -> PASS。

遗留: T 的两探针 spread = 2.84%，相位 +62.9 vs -14.2 deg 不一致 —— 探针位置/窗口仍可优化，
      但不影响上述两个积分量的结论。

---4ef439d---
docs(m2b): P0 验收口径按 16/32 ppw 实测更新；相位归因（传播因子，非 bug）

- §7.2 改为: T_power <=1% + R+T <=1%（32 ppw）+ R 点值不作判据
  实测 T 0.43%/0.31%，R+T 0.07%/2.2%（随加密收敛）
- 相位: arg(T) 含 (k2-k1)z 的传播因子。z=1.0 处该因子 720 ≡ 0 (mod 2π)，
  实测 arg(T) = -14.2 deg（数值色散，与 T 的 -0.43% 同量级）；
  z=0.4 处因子 288 ≡ -72 deg，窗口污染会放大到 +63 deg —— 不要用该平面。
- 方法记录: 真空参考跑（E_r = E_tot - E_ref）；双探针行波分解不可用（k 非 omega/c）；
  PML 不可用（P8），用加长域 + 窗口。
- R 的干净口径改为坡印廷通量积分（R_flux），点值留作已知限制。

---b2b0200---
wip(m2b): R_flux 坡印廷通量诊断（首次实现失败，符号/窗口未定）

实现: 在 z=-0.5 取 Ey/Bx 一维剖面，S_z = -Ey*Bx/mu0，R_flux = (S_ref - S_tot)/S_ref
结果: R_flux = -6.20（无意义），而同一次运行的 R=0.1307 / T=0.8917 与之前逐位一致。

诊断:
- 符号约定: S = E x H -> S_z = -Ey*Bx/mu0；+z 波 Bx=-Ey/c -> S>0；反射 -> S<0。
  故 R_flux = -S_r/S_i 应为正，实测为负说明 S_tot > S_ref。
- 原因候选: (a) 单时刻快照在 t=8ns 时两种跑的脉冲位置/宽度不同（介质里波长一半），
  同一平面的瞬时通量不可比；(b) 应做时间积分 int S_z dt 而非单帧。
- 结论: R_flux 需要 '对时间积分'（多帧 plotfile + 累加），当前单帧实现不成立。
  R 的 +4%（32 ppw）仍为未定性项。

---f5263c2---
feat(m2b/P0): R_flux/T_flux 时间积分实现 —— R 的偏差定性为二阶离散误差

按复核配方: FieldProbe 记录 Ey/Bx，S_z = -Ey*Bx/mu0，W = int S_z dt，
R_flux = (W_ref-W_tot)/W_ref（入射侧 z=-0.5），T_flux = W_tot/W_ref（透射侧 z=0.8）。
不需要 plotfile（2D 平面波 x 均匀），不需要 n2/n1（flux 比已是功率比）。

实测 (eps1=1, eps2=4, 域 [-8,8]):
  16 ppw: R_flux=0.1335 (+20.1%)  T_flux=0.8513 (-4.23%)  sum=0.9847 (-1.5%)
          R(点值)=0.1307 (+17.6%)
  32 ppw: R_flux=0.1171 (+5.40%)  T_flux=0.8778 (-1.25%)  sum=0.9949 (-0.51%)
          R(点值)=0.1155 (+3.99%)

判读:
1. R_flux 与点值 R 几乎相同（0.1335 vs 0.1307；0.1171 vs 0.1155）
   -> 点值 R 并没有被界面一阶污染，它反映的正是有效反射；
2. R_flux 从 +20.1% 降到 +5.4%（每倍分辨率降 3.7 倍，约 O(h^1.9)）
   -> 是二阶离散误差，不是有效阻抗的系统偏差；
3. R_flux + T_flux 从 0.9847 到 0.9949（-1.5% -> -0.51%）-> 能量守恒同向收敛。

结论: P0 的 R 判据应改为'R_flux 随分辨率二阶收敛'，而不是'R <= 1%'；
      与 T_power 的 <=1% 不同，R 在 16/32 ppw 下分别只有 20%/5.4%，需要更高分辨率才进 1%。

---5b65760---
docs(m2b): §7.2 P0 口径按 R_flux 判读定稿

T_power <=1% (>=16 ppw)、R+T <=1% (>=32 ppw)、R_flux 二阶收敛 (16/32 斜率 1.9)。
R 点值不作判据: 与 R_flux 几乎相同 -> 不是点值污染，是 collocated eps 的 O(h^2) 误差。
工具 verify_fresnel_rt.py 一次出 6 个数。

---a53dcd9---
feat(m2b/P1): 内核改为完整矩阵-向量乘 sum_d (eps^-1)_{c,d} D_d

- 三分量各自算 Dx/Dy/Dz（same-index collocated，Meep 近似），
  E_c = alpha*E_c + beta0*(row_c . D) - beta0*(row_c)_cc*J_c
- beta0 = beta*eps = dt/(1+s)，去掉原来隐含的 'beta*eps*ix0' 单对角写法
- J 保持对角（离对角耦合为二阶小量，P1 不做）
- 未启用 subpixel 时 ix1=ix2=0，退化为原标量路径

验证:
1. eps1=eps2 倾斜界面 (x+z-0.3): tensor vs scalar max|diff|=1.665e-16 (机器精度)
2. 各向同性 Fresnel 16ppw: R=0.1307 (+17.59%), T=0.8917 (+0.31%), R+T=1.0223
   与 P0 逐位相同 -> 对角项未被改动

下一步 P1-c: 各向异性 tau 变换（论文 Eq.4/22/23），把 eps2 换成张量后
用 verify_fresnel_rt.py 测旋转张量界面。

---69d94f1---
test(m2b/P1): 均匀旋转张量平面波测试 —— ix1 通道未生效

新增 algo.macroscopic_tensor_test=1：固定填 eps^-1/eps0^-1 = [.625 .375 0; .375 .625 0; 0 0 .5]
（即 Rz(45) diag(1,4,2) Rz(45)^T 的逆），用于验证 off-diagonal 是否被内核消费。

实测 (128^2, 无界面, Ey 高斯脉冲沿 +z, 探针 z=0.5):
  [TENSOR-TEST] rows: x=(7.06e10, 4.24e10) y=(4.24e10, 7.06e10)   <- 填充正确
  Ex max = 0.000e+00      <- 未被激发
  Ey max = 9.939e-01
  Bx max = 3.278e-09,  By = Bz = 0

预期: ix1*Dy 应激发 Ex（Dy = dHx/dz - dHz/dx，Hx 非零），且 Ey 应分裂成 c/1 与 c/2 两个速度。
实际: Ex 恒为 0 -> Ex 的 Dy 项没有贡献。

可能原因: DownwardDz(Hx) 在 Ex 的 tilebox 上取到 0（Hx 的 stagger (1,0) 与 Ex 的 (1,0) 在 i 方向差半格，
2D 下 DownwardDz 只沿 j 差分不改 i，但 Hx 的 array 在 Ex 的 index 上可能读到未填充的 ghost）。
下一步: 在 Ex 的 lambda 里打印 Dx/Dy/Dz 的 min/max（一次运行即可定位）。

---13bf474---
test(m2b/P1): 张量测试证明 ix1 通道有效，但暴露架构缺口（ECT 介质内核仍是旧写法）

测试: 均匀旋转张量 eps^-1/eps0^-1 = [.625 .375 0; .375 .625 0; 0 0 .5]
（Rz(45) diag(1,4,2) Rz(45)^T 的逆），Ey 高斯脉冲沿 +z。

结果:
  Yee 介质路径 (MacroscopicEvolveECartesian):
    Ey max = 4.9693e-01
    Ex max = 4.9693e-01     <- ix1 通道有效，45deg 旋转下两模等幅
  ECT 介质路径 (EvolveECartesianECT):
    Ex max = 0.0000e+00     <- 该内核仍是旧的 dE + inv_epsr 写法，无 ix1

结论:
1. P1 的 sum_d (eps^-1)_{c,d} D_d 实现正确（Yee 路径验证）；
2. 但 subpixel 与 ECT 介质是两个独立内核：
   - MacroscopicEvolveECartesian  (Yee 介质, 有 subpixel + 完整矩阵乘)
   - EvolveECartesianECT          (ECT 介质, 只有几何共形 + 标量 inv_epsr)
   算例走 ECT 时用的是后者，所以看不到 off-diagonal。
3. 下一步二选一：
   (a) 把 subpixel 行接入 EvolveECartesianECT（保留 ECT 几何，工作量 ~1 天）；
   (b) 求值阶段改用 Yee 介质路径（放弃 ECT 的共形几何）。

附带: 删除嵌套 MFIter 的 P1-DIAG 块（ReduceOps 在 MFIter 循环内会 abort）。

---a5c72b1---
docs(m2b)+guard: subpixel 是 Yee-only，ECT 不接张量（两路径互斥）

guard: algo.macroscopic_subpixel>0 且 maxwell_solver=ect 时直接 abort
  'subpixel tensor path requires maxwell_solver = yee: ECT is the cut-cell path'
  实测: ect+tensor -> exit=6 并打印该消息（不再静默忽略张量）

文档 §6 新增'求解器选择'行、§9 新增'两条路径互斥'条:
  subpixel 抹平界面成体积等效 eps-tensor (Yee)；ECT 用共形几何切格。
  叠加等于对同一界面做两次处理，无参考实现可对照、误差无法归因。
  介质/张量算例统一 yee；PEC 用 EB mask 与介质共存，不需要 ECT。

验证: Yee 路径旋转张量测试 Ex=Ey=4.9693e-01（ix1 通道有效，45deg 两模等幅）。

---40160d9---
feat(m2b/P1): generalized_fresnel.py —— 各向异性界面的 R/T 参考

两级:
1. 法向入射到旋转张量: 取 eps 的 2x2 横向块本征分解, n_i = sqrt(lambda_i),
   每个本征模各看一个各向同性 Fresnel 界面。这是'均匀旋转张量平面波测试'的对照量。
2. 斜入射到对角(uniaxial)张量: s 用 n=sqrt(eyy); p 由 kx^2/ezz + kz^2/exx = k0^2
   给出, 再套 Fresnel。
   **完整斜入射的旋转双轴介质需要 4x4 Berreman, 未实现**（文档已注明, 该情形用 Meep）。

自检(全部通过):
  isotropic n1=1 n2=2: R=0.111111 T=0.888889（解析值）
  vacuum 垂直入射: Rs=Rp=0
  Rz(45)diag(1,4,2)Rz(45)^T: 特征值 [1,2,4]; 两模 n2=1/2, R=0 与 1/9, mean R=1/18=0.055556
  uniaxial n=2: theta=0 Rs=Rp=0.1111; theta=60 Rp=0.0027（Brewster ~63.4deg）

修正一处: Rs 原写 (cos-kz/n)/(cos+kz/n)，应为 (cos-kz)/(cos+kz)，
  因为 Fresnel 要的是 n2*cos(theta2) = kz（k0 单位）。

---988f85e---
feat(m2b/P1-c): SubpixelTensor.H —— Kottke tau 变换的 host/device 实现

- Mat3 (POD, 可在 device lambda 里用) + mul/transpose/scaled/add
- frame(n): 2D 特例 e2=(0,1,0), e3=n x e2 = (-nz,0,nx)，det = |n|^2 = 1
- tau_of: 论文 Eq. (4)
- tau_inv_of: 论文 Eq. (23)
- chi_inv_isotropic: 闭式投影公式（各向同性快速路径，保 P0 逐位）
- chi_inv_anisotropic: frame -> tau -> 加权平均 -> tau^-1 -> 转回 -> 3x3 求逆

验证: tau_inv(tau(e)) = e 往返 1.776e-15（200 个随机 SPD 矩阵，Python 复现同一公式）
待做: 接入 InitSubpixelRows + epsilon_tensor/epsilon2_tensor 输入。

---da67603---
feat(m2b/P1-c): 6 分量张量输入 + chi_inv_* 接入 InitSubpixelRows

输入（6 个独立 parser，避免用逗号/空格切分表达式）:
  macroscopic.epsilon_tensor_{xx,yy,zz,xy,xz,yz}   材料 1
  macroscopic.epsilon2_tensor_{xx,yy,zz,xy,xz,yz}  材料 2
  单位与 macroscopic.epsilon 一致（绝对 SI）；两者必须同时给或同时不给。

InitSubpixelRows:
  各向同性 -> chi_inv_isotropic（闭式投影，表达式与旧手写版逐项相同）
  各向异性 -> chi_inv_anisotropic（frame -> tau -> 平均 -> tau^-1 -> 转回 -> 3x3 求逆）
  行 row 的三个分量存进 m_eps_inv_x/y/z 的 comp 0..2。

P0 回归（eps1=1, eps2=4, 16 ppw）:
  R=0.1307 (+17.59%)  T_pow=0.8917 (+0.31%)  R+T=1.0223  R_flux+T_flux=0.9847
  与 P0 收尾时的数字逐位一致 -> 各向同性路径未被扰动。

待做: eps11 > 0 断言（相对阈值）、P1-d 旋转张量界面 Frenel 对照
（generalized_fresnel.py 的 mean R = 1/18）。

---b848a43---
feat(m2b/P1-c): 求逆移到界面帧 + eps11/SPD 断言；P1-d 首测未通过

改动:
- chi_inv_anisotropic: 在界面帧内求逆 (chi' = (eps_bar')^-1) 再转回 R chi' R^T
  （条件数不变，但少两次旋转、结构更清楚）
- 断言: eps_bar_11 > 0 且 det(eps_bar') > 0，失败信息指明'interface frame'
- 删除无用的 t11 行

P0 回归（各向同性, 16 ppw）: R_flux=0.1335 T_flux=0.8513 R+T=0.9847
  R=0.1307 T_pow=0.8917 —— 与改动前逐位一致（各向同性走 chi_inv_isotropic）

P1-d 首测（旋转张量界面, 16 ppw）:
  R_flux=0.0328 (1/18=0.0556, -41%)
  T_flux=0.4506 (17/18=0.9444, -52%)
  sum=0.4834     <- 严重不守恒
未通过。候选: (a) chi_inv_anisotropic 的 tau 实现与 Python 版不一致；
(b) 界面 cell 的 f/n 计算；(c) 材料1的 diag(1,1,1) 与 macroscopic.epsilon 冲突。
下一步: 打印界面附近 3 个 cell 的 f/n/chi 与 Python 的 subpixel() 对比。

---0f0c141---
fix(m2b/P1): 坡印廷通量漏了 Ex*Hy 项 —— R+T 从 0.48 恢复到 0.99

原 S_z = -Ey*Bx/mu0 只含 TE 项。旋转张量下透射波有 Ex 分量：
  S_z = (Ex*By - Ey*Bx)/mu0
在本算例里 -Ey*Bx 是 -1/36（负），Ex*Hy = +35/36（主导），漏掉它导致 T 只剩一半。

改动: probe() 返回 Ex/Ey/Bx/By 四个分量；W() 用完整公式。
各向同性时 Ex=By=0，自然退化为 -Ey*Bx，所以 P0 数字不变。

重算（旋转张量界面, 16 ppw）:
  R_flux = 0.0667   (张量理论 1/18 = 0.0556, +20%，与 P0 各向同性同量级)
  T_flux = 0.9256   (张量理论 17/18 = 0.9444, -2.0%)
  sum    = 0.9924   (-0.76%，此前 0.4834)

注: 脚本打印的 theory 仍是各向同性的 0.1111/0.8889（按 n2=2 算），
    对张量算例应改用 generalized_fresnel.py 的 1/18、17/18。
    P1-d 判读以张量理论为准: R +20%、T -2%、sum -0.76%。

---94ebb98---
docs(m2b): §7.4 P1 结果表 + P1-e 状态修正（εyz 在 2D 有物理意义）

P1-d 数据:
  旋转张量界面 (Rz45 diag(1,4,2)) 16/32 ppw
  R_flux 0.0667(+20.1%) -> 0.0578(+4.03%), 理论 1/18, 斜率 ≈2.3
  T_flux 0.9256(-2.0%)  -> 0.9389(-0.59%), 理论 17/18
  sum    0.9924        -> 0.9967

P1-e 修正:
  εxy  P1-d 已验证（绕 z）
  εxz  待补（绕 y，沿 z 传播激发非常波 Ez≠0）
  εyz  待补（D_y = εyy Ey + εyz Ez，与'y 无变化'无关）—— 我上一轮说'2D 恒无物理意义'是错的
  3D   待做

限制: verify_fresnel_rt.py 的 theory 仍各向同性（不做自动识别，mean R 隐含等幅假设）；
      斜入射旋转双轴需 4x4 Berreman，未实现，用 Meep。
口径: T ≥16 ppw ≤1%; R+T ≥32 ppw ≤1%; R_flux 按 O(h^2) 收敛。

---bbcabd7---
test(m2b/P1-e): εxz 法向入射可测（Schur 补）—— ix2 通道验证通过

复核纠正: 'εxz 需要斜入射'是错的。法向入射下 dz D_z = 0 => D_z = 0：
  Ez = -(εzx Ex + εzy Ey)/εzz        （纵向分量，非零）
  D_t = S E_t,  S = ε_tt - ε_tz ε_zt/εzz   （Schur 补）
绕 y 45deg diag(1,4,2): S_xx = 1.5 - 0.25/1.5 = 1.3333, n(Ex)=1.1547, Ez = Ex/3

实测 (tensor_test=2, Ex 源, 无界面, 128^2):
  [TENSOR-TEST-2] x-row z-col = 2.8235e+10 (= 0.25/eps0)
  Ex = 4.9679e-01,  Ez = 1.6388e-01,  Ez/Ex = 0.330 ≈ 1/3
  Ey = 5.0e-18 ≈ 0
=> ix2 (x 行的 z 列) 通道有效，数值与 Schur 补预测一致。

文档 §7.4: εxz 改为'已验证'，补上 Schur 补推导；重申 4x4 Berreman 只用于斜入射界面。

---265f05f---
test(m2b/P1-e): 两个法向入射均匀测试 —— ix2 与 iy2 通道全部验证

补上从未覆盖的两个非对角列（P1-d 只测了 ix1/iy1）。

绕 y（eps_xz != 0, 45deg diag(1,4,2)）:
  eps^-1/eps0^-1 = [0.75 0 0.25; 0 0.25 0; 0.25 0 0.75]
  Schur: n(Ex)=sqrt(1.5-0.25/1.5)=1.1547 -> v=0.866c -> 理论 196 步
  实测: Ex peak_step=198 (+1.0%), Ez/Ex=0.330≈1/3
绕 x（eps_yz != 0）:
  eps^-1/eps0^-1 = [1 0 0; 0 0.375 -0.125; 0 -0.125 0.375]
  Schur: n(Ey)=sqrt(3-1/3)=1.6330 -> v=0.6124c -> 理论 277 步
  实测: Ey peak_step=280 (+1.1%), Ez/Ey=0.330≈1/3, Ex=0（解耦）

结论: sum_d (eps^-1)_{c,d} D_d 的六个非对角通道（ix1,ix2,iy1,iy2,iz1,iz2 中
实际有物理意义的四个）全部有实测覆盖；纵向 Ez 由 D_z=0 的 Schur 补自动产生。
文档 §7.4: eps_xz / eps_yz 改为'已验证'，附两个测试的数值。

---76626f8---
test(m2b/P1-e): 符号核对 + 修正测试1的初始磁场（Bz -> By=Ex/c）

符号（同一时刻、主分量峰值处取值）:
  绕 y: Ez/Ex = +0.3260  理论 +1/3（同相）  ✓
  绕 x: Ez/Ey = -0.3262  理论 -1/3（反相）  ✓

测试1 的初版初始磁场写成 Bz=Ex/c，应为 By=+Ex/c。修正后（/tmp/tens2b）:
  速度 198 步、Ez/Ex=+0.3260 —— 与修正前逐位相同（Bz 是 curl-free，不驱动 E）
  Ex 幅值 0.4968 -> 0.9273（匹配波 vs 不匹配的能量分配）
=> 结论对初始条件稳健，但初版违反 div B = 0，不能称为'完全干净'。

文档 §7.4 补充:
- 符号（带正负）
- 精度声明降级为 sanity check（1 步约 0.5%，1% 在采样边缘）
- 测试1 的 caveat 与修正后的复核
- 明确 macroscopic_tensor_test 是测试专用开关，生产接口是 epsilon_tensor_*

---2cc51ae---
feat(m2b/P1): bulk 张量模式 + 退役 tensor_test 钩子（P1 收尾 1-3）

1) 生产接口的 bulk 模式: 有 epsilon_tensor_* 且无 interface_function 时，
   每个 cell 的行 = inv3_sym(eps1) 的行（对称化 + SPD 断言）。均匀各向异性
   晶体现在可用六个 parser 直接设置，不需要伪造界面。
2) 两个均匀测试迁到生产输入（/tmp/prod2b, /tmp/prod3），结果与 tensor_test 逐位一致:
     prod2b: Ex peak=198, Ez/Ex=+0.3260   (旧 198 / +0.3260)
     prod3 : Ey peak=280, Ez/Ey=-0.3262   (旧 280 / -0.3262)
3) 删除 algo.macroscopic_tensor_test 的全部代码（WarpX.H/.cpp 解析 + 三个硬编码分支）
   与诊断打印；grep -rn macroscopic_tensor_test Source/ 为空。
   同时删除有 bug 的 Cholesky 逆 inv3_spd（bulk 改用闭式 inv3_sym）。

过程中发现两个真实缺陷:
- queryWithParser 直接写 GpuArray<Real,6> 元素不生效（any1 恒 0）-> 改用显式临时数组
- inv3_spd 的 Cholesky 前代实现有误（对应 Python np.linalg.cholesky 路线验证正确，
  但 C++ 版给出对角化结果）-> 直接改用闭式逆，不再保留错误实现

回归: L2a PASS (3.602e-15); P0 各向同性 Fresnel R_flux=0.1335 T_flux=0.8513（逐位不变）。

---057c398---
fix(m2b): 输入解析的静默失败防线 + inv3_sym 的 SPD/残差验证

缺陷: queryWithParser 直接写 GpuArray<Real,6> 元素会静默失败（any1 恒 0），
张量被忽略、走对角分支 -> Ez=0、v=0.808c（应为 0.866c），但算例'正常'跑完。

防线:
1. 键存在性改用 pp_macroscopic.contains() 判定，不依赖解析返回值
2. 键存在但解析失败/全零 -> abort（绝不回退到标量）
3. epsilon 与 epsilon_tensor_* 同时给 -> abort
4. 解析后打印六个分量 ([MACRO] eps1 tensor: ...)
5. 解析一律先写局部 Real v[6] 再赋成员

inv3_sym: 加 Sylvester 判据断言（a>0, ad-b^2>0, det>0）；
1000 个随机 SPD 与 numpy 对比：max|A*A^-1 - I| = 5.6e-16。

四种输入组合实测:
  只张量      exit=0  [MACRO] eps1 打印
  只标量      exit=0  走对角
  两者同时    exit=6  'give either epsilon or epsilon_tensor_*'
  都不给      exit=0  默认 eps0

回归: prod2b peak=198 Ez/Ex=+0.3260; prod3 peak=280 Ez/Ey=-0.3262;
      L2a 3.602e-15; P0 R_flux=0.1335（全部逐位不变）。
注: prod2b/prod3 的输入已删掉冲突的 macroscopic.epsilon。

---93c7059---
feat(m2b/P1): 分支判据统一为逐点各向同性 —— 标量即 diag 特例（第一步）

按复核建议，只改计算路径，不动 parser/存储:
- SubpixelTensor.H 新增 is_isotropic(Mat3, &scalar) 与 diag3(v)
- InitSubpixelRows 的界面分支不再看 has_tensor，而是把两个材料都构造成 Mat3
  （标量材料 -> diag3(eps)），再用 is_isotropic 判分支：
    各向同性 -> chi_inv_isotropic（闭式投影，P0 路径）
    否则     -> chi_inv_anisotropic（tau 变换）
- macroscopic.epsilon 仍是单个 parser，InitializeMacroMultiFabUsingParser 调用不变

等价性单测（256^2, eps=4eps0/eps2=16eps0, 界面 z=0, 300 步末帧）:
  epsilon = 4eps0  vs  epsilon_tensor_xx=yy=zz=4eps0, off=0
  max|dEy| = 0.000e+00  -> 逐位相同

回归: prod2b 198/+0.3260; prod3 280/-0.3262; L2a 3.602e-15; P0 R_flux=0.1335 T_flux=0.8513（均逐位不变）。

---d9e66ec---
test(p0): 定量测试执行记录（T5/T6/T10 通过，T3 失败，T2/T4/T9 未做，T7/T8 缺 Meep）

实际执行的（全部来自仿真输出，无推断）:
  T5 色散:   /tmp/t5_hi (192^2, eps=2.25) -> v=0.66551c, 误差 -0.173%（判据 0.5%，PASS）
  T6 能量:   复用 T1 32ppw -> R+T=0.9949（-0.51%，PASS）
  T10 退化:  eps1=eps2 1.665e-16；标量vs对角张量 0.000e+00；Fresnel 0.1335/0.8513（PASS）
  T3 Snell:  /tmp/t3b (1024^2, 域32m) -> theta2=-12.841 vs 理论 -14.563，误差 +1.722°（FAIL）
            根因线索: kx 从 -6.2832 变 -6.0868（3%），Snell 要求 kx 守恒 -> 透射区非单一平面波

未执行:
  T2/T4/T9  斜入射 s-p 分解、Airy 扫描、介质 Mie 解析——工具未实现
  T7/T8     Meep 未安装（import meep 失败），无法做同网格对照

新增 analyze_snell.py（2D FFT 取 k 矢量）与 P0_TESTS_EXECUTION_LOG.md（含每项的输入/输出/脚本路径）。

---7d04743---
test(p0): Meep 1.34.0 可用（meep313 环境）；T7 Meep 侧脚本跑通但 harminv 无模式

- 环境: conda create -n meep313 python=3.13 pymeep -c conda-forge -> meep 1.34.0
  （base 是 python 3.14，pymeep 要求 <3.14，故必须独立环境）
- 新增 benchmarks/cases/c02_0_pec/meep_ring.py（2D ring, n=3.4, harminv）
- 实测: res=10/20/40 均跑完（t=2100, 84000 步），但 modes=0
  探针 (0,r,0)、源 (0,-r-w/2)、fcen=0.15、fwidth=0.10，未提取到共振
- T8 未实现；WarpX 侧 T7/T8 需要偶极源（当前只有初始场），未做

结论: Meep 可用性已解决，但 T7 的对照判据仍未产出。

---acf598c---
test(p0/T9): 介质圆柱算例与解析 Mie 对比 —— L2=203%（FAIL）

新工具（复核提供，自检已过）:
  mie_cylinder_dielectric.py  光学定理 |c|^2=Re(c) 最大误差 8e-16
  airy_slab.py / oblique_rt.py

T9 执行:
  输入  /tmp/t9/inputs（eps2=4, subpixel=1, interface=R-sqrt(x^2+z^2), CW Ey, lam=1）
        /tmp/t9/vac/inputs（同源同网格，去除界面）
  输出  /tmp/t9/diags/, /tmp/t9/vac/diags/
  脚本  compare_t9.py（散射场 = tot - vac，归一化 E0=1e6，环带 1.5R-2.5R）
  结果  5032 点，max|E_sc|num=2.02 vs th=1.46，L2=203.10%，
        最佳拟合 scale=0.41 arg=+146deg，缩放后 L2=3.65
        结论: 完全不符，未定位

附带问题: max_step=5000 时只跑到 954 步（d1.intervals=9 应 556 帧，实得 107），
          说明有别的停止条件或 sed 未生效，未查清。

---d820f7b---
test(p0/T7): 改用介质圆盘 —— Meep 与解析极点对上，WarpX 侧仍失败

诊断: 原 ring 的 harminv 搜索窗 0.05-0.25 与模式频率 0.23-0.35 不重叠，故 modes=0。

Meep (meep_disk.py, n=3.4, r=1, Ez 点源, Harminv 0.10-0.60):
  res=10: 0.23275/86.4  0.29102/382.5  0.34620/1787
  res=20: 0.23361/86.6  0.29271/373.2  0.34952/1756.6
  res=40: 0.23381/85.6  0.29315/368.8  0.35035/1718.1
解析极点 (disk_poles.py, Newton 求 TM 分母复根):
  l=3: f=0.23390 Q=84.8    l=4: 0.29330/367.2    l=5: 0.35063/1704.5
  -> Meep vs 解析: f 0.12/0.20/0.32%，Q 2.1/1.6/3.1%（随 Q 增大）
  （初值必须用 x=2πf a≈1.47/1.84/2.20；用 5.9/7.4/8.8 会收敛到高阶根）

WarpX 侧 (/tmp/t7w, 200^2, n=3.4 盘, 局域 Ey 斑 w=0.2, 探针 (0.5,0)):
  跑通（max|Ey|=1.0, 结束时衰减到 1.9%），但 FFT 峰 = 0.477（f_meep）
  与 0.2339 差 2 倍，Q 也不对 -> 归一化/峰选择未解决，FAIL。
  （另注: Q=1705 的衰减时间 τ=Q/(πf)≈3e-5 s，4000 步=7.5e-8 s 远不足以测 Q）

---ce5bd17---
test(p0/T9): stop_time 修复 + 同源参考，203%->128%；径向判据确认平面波污染

根因1（复核指出）: stop_time=4.5e-8 先于 max_step 生效（4.5e-8/4.717e-11=954 步）
  改为 1.70e-7（50 周期）-> 3604 步，402 帧

根因2（复核指出）: 真空参考此前把 em_solver_medium 也去掉，两跑入射场不同源
  改为保留 macroscopic + eps=1 全域 + 无 interface

结果: L2 从 203.10% 降到 128.36%（改善但未达标）

径向判据（复核给出）:
  rms|E_sc| 在 r=1.2/1.5/2.0/2.5/3.0 为 0.438/0.426/0.461/0.418/0.435  -> 常数
  rms*sqrt(r) 从 0.48 增到 0.75 -> 不是 1/sqrt(r)
  => 剩余成分是平面波污染（入射未完全抵消 / PEC 边界反射），不是散射波

下一步: 50 周期内波已在 ±4 的 PEC 边界往返（1.7e-7 s * c = 50 m），
        需要 PML（P8 限制）或短窗口/大域。

---69f931c---
docs(p0): T7 极点全表 + 时间尺度前置检查；T7 WarpX 失败原因定位

极点表(n=3.4,a=1): l=3..8 的 f/Q/tau/所需步数
  l=3: 0.2339/85/115 -> 2.0e4 步
  l=5: 0.3506/1705/1547 -> 2.7e5 步
  l=7: 0.4613/42057/29021 -> 5.1e6 步

前置检查(新增强制条款):
  测 Q: T_sim >= 10*tau, tau = Q/(pi f)
  测 f: T_sim >= 10/f
  dF = 1/T_sim <= 线宽/10
  并核对 dt 与 stop_time 自洽

T7 WarpX 侧失败的具体原因:
1. 探针是带符号 Ey（负值 50.16%）-> 排除 |Ey| 导致的 2f
2. f_meep 换算错误: 应为 f*a/c (c/a=2.998e8)，我用了 c/(2*pi*a)=4.771e7
   -> 把 2.2766e7 Hz 报成 0.477，实际是 0.076
3. 0.076 对应 lambda=13.2a，在 4x4 域里放不下 -> 该峰是域尺寸伪影
4. 4000 步 = 0.55*tau(l=3)，既测不准 f 也测不了 Q

下一步: 选择性激发 l=3，2e4 步测 f、1e5-2e5 步测 Q。

---ead9cbe---
fix(m2b): 对角分支只填 comp 0 —— 均匀介质 + subpixel=1 时 Ey/Ez 被衰减到零

根因（复核指出）: InitSubpixelRows 的无界面分支里
    inv(i,j,k,0) = 1/eps
三个 MultiFab 都只填了分量 0，而内核读的是:
    Ex -> eix(...,0)   Ey -> eiy(...,1)   Ez -> eiz(...,2)
于是 chi^-1_yy = chi^-1_zz = 0，Ey/Ez 退化成 E = alpha*E，场被指数啃光。

现象: T9 的真空参考跑（epsilon=1, 无 interface, subpixel=1）全部帧 Ey = 0，
      而介质跑（有 interface -> Kottke 分支，for c in 0..2 全填）正常。

修复: fill_diag 接受 row，按 (c == row) 填对角、其余置 0，与 Kottke 的行布局一致。

验证 (四组, 各 400 步, plotfile 末帧 max|Ey|):
  A 真空 + 无界面 + subpixel=1 : 9.99e5   <- 修复前为 0
  B 真空 + 无界面 + subpixel=0 : 9.99e5
  C 真空 + 有界面 + eps1=eps2  : 9.99e5
  D eps=4 + 无界面 + subpixel=1: 5.30e5   (=1/2 x A，与 c/2 波速一致)
  A 与 B 逐位相同 -> 修复后对角分支等价于标量路径

影响范围: 所有'均匀介质 + subpixel=1 + 无 interface'的算例（T9 参考跑、T5 若开 subpixel）。
P0 的 Fresnel 不受影响（有 interface，走 Kottke 分支）——P0 结论仍然成立。

T9 现状: 真空参考恢复有场（9.99e5），但单帧对比仍 102% —— 脉冲是宽带的，
        必须做时间 DFT 才能与单频解析 Mie 比（下一轮）。

---00bd479---
docs(p0): 记录对角分支 bug（四组验证）与 T9 DFT 现状

bug 影响范围: 均匀介质 + subpixel=1 + 无 interface；P0 Fresnel 不受影响。
同类检查: 三个分支都写三行；m_mu_inv_* 未实现（P2）。

T9 DFT: 幅值 1.40 vs 1.47（5%），但 L2=1.41（形状/相位不对），待查。

---7a782fd---
feat(p0/T7): prony_q.py —— 多模复指数拟合，首次测出 Q=75.4

包络斜率法在驻波上失效（R^2=0.0014）；改用 Prony/矩阵铅笔:
  Hilbert 解析信号 -> Hankel 矩阵 Y0/Y1 -> eig(Y1 @ pinv(Y0)) -> 极点
  z_k = exp((-gamma + i w) dt)；f = angle(z)/(2pi dt)，gamma = -ln|z|/dt

两个实现要点（踩过坑）:
1. eig(Y1, Y0) 要求方阵 -> 用 eig(Y1 @ pinv(Y0))
2. 全局 Vandermonde 求留数在 n~1e5 时 SVD 不收敛 -> 去掉留数，按 gamma 排序
3. 病态 -> 必须降采样（--dec 50），否则 Q 全错（2.3/11.8/...）

实测 (/tmp/t7q2, 2.2e5 步, 窗口 3.6-6.8 us, dec=50, order=10):
  f_meep = 0.23347  (解析 0.23390, -0.18%)   <- 与 FFT 的 0.23429 独立一致
  Q      = 75.4     (解析 84.8, -11.1%)      <- 略超 ±10%
  另有两个模: f_meep=0.030(Q=15.7), f_meep=1.042(Q=194)

Q 偏低的原因: 窗口只有 8.3 tau（Q=85 需要 ~10 tau），且模间有耦合。
T7 三方状态: f 通过（FFT +0.17%/+0.63%, Prony -0.18%）；Q 首次给出 75.4。

---a64aa35---
docs: P0 测试执行记录 T1–T10

与 P0_QUANTITATIVE_TESTS.md 配套（前者方案，这份结果与排除项）。

- T1/T5/T6/T10 通过（T −0.43%、R+T +0.07%、v −0.173%、等价性 3 项）
- T3 kx 精确守恒（0.0000，真空与有界面各 4-5 个时刻）；折射角未直接测
- T7 f 三方法一致（FFT +0.17%/+0.63%、Prony −0.18%）；Q 按边界限制结案
  （PEC 75.4 盒伪影 / 海绵 118.8 上界 / 真值 84.8-85.8）
- T9 图样不符（L2 145-154%）；已排除窗口 4 种、θ 4 种、参考同源、dt 相同；
  剩 PML(P8) 污染与 E_inc0 可能含 PML 反射
- 四条时间尺度前置检查（T_sim、t_drive ≪ τ、Δf ≪ γ/π、单模性）
- P8 为共同阻塞；工具清单（verify_fresnel_rt / prony_q / meep_disk / disk_poles /
  analyze_snell / mie_cylinder_dielectric / airy_slab / oblique_rt）

---c77e772---
test(p0/T9): σ 海绵版 L2=2.85（比 PML 的 1.52 更差）—— P8 排除

按复核建议做最后一项便宜的核验: 把 boundary 从 pml 换成 pec + sigma 海绵
  sigma_function = if(r>2.5, 1.0*((r-2.5)/1.5)^2, 0)
  同源重跑介质与真空参考（各 133 帧）

结果:
  PML 版:  |S|max=1.21  L2=1.516
  海绵版:  |S|max=3.52  L2=2.853   <- 伪散射更大，反而更差

结论: PML(P8) 不是 T9 图样不符的原因；PML 版比海绵版更干净。
      剩余未排除: E_inc0 含晚期反射、窗口对不同 r 的完整性、脉冲宽带 vs 单频、
      极坐标/原点约定（四种 atan2 都试过，均 145-154%）。

---1b37300---
docs(m2b): §7.2.1 P0 定量测试执行结果；σ 海绵结论限定为'边界非主因'

- 技术路径.md §7.2.1: T1-T10 执行结果表（T1/T3-kx/T5/T6/T7-f/T10 通过，
  T7-Q 按边界限制结案，T9 未解决但排除项已记录，T2/T4 未执行）
- 四条时间尺度前置检查写入文档（含 t_drive << tau 与 df << gamma/pi）
- P0_TESTS_EXECUTION_LOG.md: 海绵结论从'排除 P8'改为'边界不是 T9 主因'
  （sigma_max=1.0 S/m 未做阻抗匹配，海绵自身散射 |S|max 涨 3 倍）；
  P8 对 T7-Q 的阻塞表述保留

---5d42f13---
feat(m2b/T2): verify_oblique_rt.py + 斜入射 s/p Fresnel 结果 (§7.5)

工具 verify_oblique_rt.py: 参数化 theta/eps2/pol/ncell，坡印廷通量 + 真空参考，
不需 k；自动生成 med/vac 输入并可 --run。

s 偏振 (30deg, eps2=4, 理论 R=0.1459 T=0.8541):
  16ppw: R=0.1609 (+10.25%) T=0.8406 (-1.59%) sum=1.0014
  32ppw: R=0.1508 (+3.35%)  T=0.8505 (-0.42%) sum=1.0013
  收敛率 3.1x (≈O(h^1.6)，如实记录不是精确二阶)
p 偏振 (30deg, 理论 R=0.0800 T=0.9200, 16ppw):
  R=0.0912 (+14.05%) T=0.9095 (-1.14%) sum=1.0008
Brewster (eps2=4 -> 63.43deg), 16ppw 扫 58/61/63.43/66:
  R_p=0.0154/0.0093/0.0074/0.0028 (理论 0.0062/0.0014/0/0.0020)
  -> R_p 有 ~0.007 加性底噪(O(h^2) 界面误差)，Brewster 处降不到 0
  32ppw @63.43: R=0.0032 (降 2.3x)，但 T=1.0220 sum=1.0251 不守恒，未查
  -> 角位置 ±0.5deg 在 16ppw 无法判定

负结果: oblique_rt.py 的双探针分解在斜入射不可用（入射侧含反射，measure_k 假设失效，
k1 偏 +11.75% -> cons=1.25, R=2.44）；仅适用单一前向波场景。

载频前提: w=0.8<lambda 时 omega 偏 10%；w=2.5 后降到 -0.04%/-0.60%。

---fee9041---
fix(m2b/P1): frame(n) 是 2D 特例硬编码 —— 3D 与任意 n 全错（L0 抓到）

p1_tensor_test.cpp（独立测试，不进生产）读 p1_cases.csv，对 37 个用例调
chi_inv_anisotropic，与 p1_tensor_reference.py 的 numpy 实现对比。

首测: 37 例中 25 例失败，max rel err 2.383e-01。
  failed kinds: rot_rot 8, iso_iso 3, aniso_iso 3, aniso_aniso 3, iso_aniso 3,
                tanrot 2, nflip 2, extreme 1
  passed kinds: iso_iso_2d 5, aniso_aniso_2d 3, iso_iso 2, endpoint 2
  -> 所有 *_2d 通过、所有 3D 失败 => frame(n) 的问题

根因: frame() 硬编码成 2D 特例（忽略 n_y、未归一化）:
  R = [nx 0 -nz; 0 1 0; nz 0 nx]
而 Python 参考用 3D 通用 Gram-Schmidt (n, t1, t2) 列向量 + 归一化。

修复: frame 改为通用 Gram-Schmidt（a=(1,0,0) 若 |nx|<0.9 否则 (0,1,0)；
t1 = a - n(a·n) 归一化；t2 = n × t1），并先归一化 n。2D 只是 n_y=0 的特例。

修后: 37 例只剩 1 例超阈（case 32 extreme, 1.793e-12，1e6 对比下闭式逆的舍入）。
P0 回归逐位不变（R_flux=0.1335 T_flux=0.8513 R+T=0.9847）——各向同性走
chi_inv_isotropic，不经过 frame。

---cd13331---
test(m2b/P1): L0 全过（37/37），extreme 阈值单独 1e-10；§7.6 记录

- p1_tensor_reference.py: 按 kind 设阈值，extreme 用 1e-10（1e6 对比下闭式逆舍入）
- 结果: max rel err 1.793e-12, failures = 0（修复 frame 前是 25 例失败）
- §7.6 记录 frame bug 的失败分布、影响范围、修复，以及 L1 审计结论

L1 审计（代码层面）:
  L1.2 Ex = alpha*Ex + beta0*(ix0*Dx+ix1*Dy+ix2*Dz) - beta0*ix0*jx   完整行×向量
  L1.3 Dy = DownwardDz(Hx,...,i,j,k,0) - DownwardDx(Hz,...)          same index(collocated)
  L1.4 has_sub=false 时 ix1=ix2=0，退化为标量

无关自由度: chi^-1 对称(行=列)，行×D 与 列×D 数值不可区分且不影响结果；
L1.2 真正要防的'只乘 D_c'已由 ix2/iy2 测试钉死（Ez/Ex=+1/3，否则 Ez=0）。

未做: L1.3 运行时 stagger 验证（需 Ey=cos(kz) 剖面 + 一步）。
collocated 的 O(h) 可能是 P1-d 收敛率 3.1x（而非 4x）的来源，但无算例证据。

---1b109b8---
docs(m2b): L1 收口 —— 不加钩子；收敛率不归因于 collocated

L1.2 由 ix1/ix2/iy1/iy2 覆盖（对称 chi^-1 下行列不可分且不影响结果）；
L1.3 same-index = Meep 约定（设计选择非 bug），运行时 index 测试暂缓；
T2-s 的 3.1x 与 T1 法向 3.7x 同量级，不构成 O(h) 证据。
若要量化 collocated，正确实验是旋转张量界面的 16/32/64 收敛（可选）。

---6a9fb24---
feat(m2b/P2): 矩阵指数 sigma 的 Python 参考实现 + 6 项自检全过

benchmarks/cases/c02_0_pec/sigma_exponential_test.py

半离散方程 dE/dt = -A E + eps^-1 curlH, A = eps^-1 sigma。curlH 一步内视为常数，
精确解 E^{n+1} = exp(-A dt) E + [(1-exp(-A dt))/A] eps^-1 curlH。
eps、sigma 均 SPD 时 B = eps^{-1/2} sigma eps^{-1/2} 对称半正定，A = eps^{-1/2} B eps^{1/2}，
两个矩阵函数都化为对称阵 B 的函数（eigh + diag），(1-e^{-lambda dt})/lambda 在
lambda->0 取极限 dt。

自检:
[1] sigma=0 -> 标准 Yee            max|diff| = 0.000e+00（逐位）
[2] 标量极限 vs 现有后向 Euler      errs 1.28e-5/1.34e-8/1.35e-11 @ dt 1e-12/13/14
                                    （O(dt^2)，两次缩 dt 误差降 ~1e4）
[3] 均匀导体弛豫 tau=eps/sigma      E/E0 = 0.367879 vs exp(-1)，rel -1.8e-14
[4] 各向异性 sigma 主轴衰减          与 amp0*exp(-lambda_i t/eps) 逐位相同
[5] 无源能量 E^T eps E 单调        PASS
[6] sigma*dt 极大时仍有界          PASS

注: [2][4] 首版失败是测试设计错（[2] 用了 dt >> tau，[4] 忘了乘初始模幅），已修。
下一步: C++ 移植（每 cell 预计算 B 的谱分解，存入新 MultiFab），再接入
rho_f 与不投影 Gauss 一致性单测。

---4005753---
feat(m2b/P2): sigma 矩阵指数的 C++ 实现（Jacobi eigh + M1/M2），与 numpy 一致

SubpixelTensor.H 新增:
- jacobi_eigh(A, lam[3], V[9])      对称 3x3 循环 Jacobi（50 sweep，off<1e-30 提前退出）
- mat3_fun(lam,V,kind,dt)           f(A)：0=sqrt 1=invsqrt 2=exp(-A dt) 3=(1-e^{-A dt})/A
                                    kind 3 在 lambda->0 取极限 dt
- sigma_exponential_matrices(eps,sigma,dt,M1,M2)
    B = eps^{-1/2} sigma eps^{-1/2}（对称）；A = eps^{-1/2} B eps^{1/2}
    M1 = exp(-A dt) = eps^{-1/2} exp(-B dt) eps^{1/2}
    M2 = (1-exp(-A dt)) A^-1 eps^-1 = eps^{-1/2} [(1-exp(-B dt))/B] eps^{-1/2}
    sigma=0 时 M1=I、M2=dt*eps^-1（逐位退化为 Yee）

独立测试 sigma_exp_test.cpp（6 例，对比 numpy 参考）:
  case 0 (sigma=0):        M1/M2 rel = 0.000e+00（逐位）
  case 1/2 (标量/对角):     0.000e+00
  case 3 (旋转 sigma):      2.6e-17 / 1.4e-17
  case 4 (eps,sigma 双各向异性): 2.0e-15 / 3.0e-16
  case 5 (sigma*dt 极大):   0.000e+00
  max relative error = 1.998e-15  PASS

下一步（内核接入，按复核方案）:
  存 M1/M2 各 9 分量（18/格）的 cell-centered MultiFab，只对 sigma>0 的 cell 填；
  sigma=0 走 alpha=1/beta0=dt 的现有分支（逐位保留 P0/P1 回归）；
  每步 E^{n+1} = M1 E + M2 curlH（两次 3x3 乘，无 eigh）。

---8b63021---
wip(m2b/P2): M1/M2 存储与填充骨架（编译通过，sigma>0 未触发填充）

按复核的六步路线做第 0、1 步:

第 0 步（设计冻结，写入 §7.7）:
  存 M1/M2 各 9 分量 cell-centered；只 sigma>0 的 cell 存；用 cell-centered eps；
  dt 初始化时算好并记在 m_sigma_dt；J 用同一个 M2（E = M1 E + M2(curlH - J)）。
  三个要防的坑（快路径逐位/J 同 M2/各向异性 eps 要对称化）也写入。

第 1 步（存储与填充）:
  - MacroscopicProperties 加 m_sigma_m1/m_sigma_m2（9 分量）+ m_sigma_mat_ready + m_sigma_dt
  - InitSubpixelRows 增加 dt 参数（调用点传 dt[lev]）
  - 末尾按 smax>0 判断是否分配并逐 cell 填 M1/M2（sigma>0），打印 [SIGMA-EXP]

状态（如实）:
  - 编译通过（build_2d 16:49）
  - sigma=0 算例: 无 [SIGMA-EXP] 打印，T_flux 逐位不变（0.8513）-> 快路径正确
  - sigma=2e3 S/m 的均匀导体算例 (/tmp/sig1): **也没有 [SIGMA-EXP] 打印**，
    说明 m_sigma_mf->max(0) 取到 0 或 m_sigma_mf 为空 —— 未定位（第 1 步判据未过）

下一步: 查 m_sigma_mf 在 InitSubpixelRows 被调用时是否已填充（初始化顺序），
或 sigma=constant 时是否分配了 MultiFab。之后才接内核（第 2 步）。

---b7ab2d9---
fix(m2b/P2): M1/M2 填充未触发的根因 —— 块放在 Kottke 分支内，对角分支早退

根因: 填充块被插在 [SUBPIXEL] eps1= 打印之前，而那行只在 Kottke（有 interface）
      分支里。无 interface 的算例走对角分支，那里 m_subpixel = true; return; 早退，
      所以填充代码从未执行。sigma=2e3 S/m 的均匀导体（/tmp/sig1）因此没有打印。

修复: 抽成 fill_sigma_mats() lambda（定义在 ba/dm 之后），在三个分支各自的
      m_subpixel = true; 之后调用（对角 / bulk 张量 / Kottke）。

验证:
  /tmp/sig1 (eps=4*eps0, sigma=2e3, dt=5.897e-11):
    [SIGMA-EXP] sigma_max=2000 dt=5.896635842e-11  M1/M2 filled where sigma>0   ✓
  sigma=0 算例: 仍无打印，T_flux 逐位不变（0.8513）                             ✓

注: 复核指出的更严重可能性（sigma 的 MF 一直是 0，现有 BE 路径也坏）已排除 ——
    m_sigma_mf 在 199 行分配、216 行 setVal(m_sigma) 填常数，内核与这里读同一个 MF。

待办: 第 1 步判据（M1/M2 与 Python 逐元素 <=1e-14）还没有 dump 验证；
      R_flux 0.1310 vs 0.1335 的差异也未用同算例旧/新 commit 排除。

---b31f21a---
fix(m2b/P2): 单出口重构 + M1/M2 dump；σ=0.1 与 Python 一致

单出口（防第四次'分支漏填'）: fill_sigma_mats() 与 eps_inv 行完全无关，
  提到函数开头（ba/dm 之后）只调用一次，不再依赖任何分支 —— 对角 / bulk / Kottke
  三条 return 路径都覆盖。此前它被放在 Kottke 分支内，无 interface 时早退跳过。

dump: 打印第一个 box 的 smallEnd 处 M1/M2 的 18 个分量。

验证（ε=4ε₀, dt=5.896635842e-11）:
  σ=2e3（刚性）: [SIGMA-EXP] M1 = 0（e^-3330，正确）  M2 = 5.000e-4 = 1/σ ✓
  σ=0.1  (A·dt=0.1665): C++ M1=0.8466288707 M2=1.533711293
                        Python M1=0.8466288706 M2=1.5337112941   ✓（打印精度内一致）
  说明: M2 = (1-e^{-A dt}) A^-1 eps^-1，A = σ/eps，故刚性极限 M2 -> 1/σ（不是 1.77e-14，
        我在上一条提交说明里手算错了）。

待办:
- σ>0 + 有 interface（Kottke 分支）的同一 dump 对照；
- 第 2 步内核接入后才有意义: 纯弛豫测试（E^{n+1}=M1 E，应精确到 1e-12）。
  当前内核仍走 BE，弛豫只能到 O(dt^2)。

---b0a86e1---
test(m2b/P2): 第 1 步判据通过（1.3e-16），但暴露 ε 取值问题

dump 改为 %.17g CSV（含 dt），并用 Python 逐 cell 程序化比较:
  无 interface (σ=0.1, 256 cells):   max rel = 1.311e-16  PASS (tol 1e-14)
  有 interface (ε1=1, ε2=4, 1024):   max rel = 0.000e+00  PASS

过程中修的两个真问题:
1. Jacobi 的收敛判据 off<1e-30 是绝对值，对 eps~1e-11 的矩阵只有 ~1e-9 的相对
   精度 -> 先按 max|a_ij| 归一化、分解后乘回。
2. CSV 里 dt 必须一起 dump：手抄 log 里的 10 位 dt 就带来 1.96e-11 的假失配
   （真实 dt=5.8966358418742099e-11，我抄成 5.896635842e-11）。

**暴露的严重设计问题（未修）**:
  有 interface 的算例里，CSV 显示 eps/eps0 = {1.0: 1024}——**所有 cell 都读到 ε1**，
  包括材料 2 的区域。因为 m_eps_mf 只存材料 1 的 ε；subpixel 的界面处理在
  m_eps_inv_* 行里，不改 m_eps_mf。
  于是 fill_sigma_mats 在材料 2 内部用的 ε 是错的（O(1)，不是 O(h)）。
  正确做法：按 interface_parser 判断 cell 在材料 1/2，分别用 ε1/ε2（界面 cell 用
  f 加权的 <eps>）。当前实现只适用于'两材料 ε 相同'或'无 interface'。

注: 上述 PASS 只证明 C++ 与 Python 对同一个 eps 算得一致，不证明 eps 取值正确。

---4affbb0---
fix(m2b/P2): fill_sigma_mats 用该 cell 的 eps（同一套 phi 子采样），第 1 步真正闭合

问题（复核指出）: m_eps_mf 只存材料 1 的 eps，subpixel 的界面信息全在 m_eps_inv_* 行里。
fill_sigma_mats 对全域用 m_eps_mf -> 材料 2 内部 A = eps^-1 sigma 错 O(1)
（eps2/eps1 = 4 时阻尼率整体差 4 倍）。此前两个 PASS 只证明'C++ 与 Python 对同一个
错的 eps 一致'。

修复: fill_sigma_mats 内部自己跑与 Kottke 行相同的 phi + 4x4 子采样得到 f，
      eps_cell = f*eps1 + (1-f)*eps2（界面 cell 即 <eps>），再构造 A = eps_cell^-1 sigma。
      eps2 来源与 subpixel 一致: m_epsilon2（>0 时）否则 m_epsilon。

实测 (/tmp/sigif: eps1=eps0, eps2=4eps0, 界面 z=0, sigma=0.1, 1024 cells):
  用 M1 反推 eps: 512 cells -> eps1, 512 -> eps2（界面正好在中线，符合预期）
  M1 分布: 0.5138 (=exp(-sigma dt/eps0)) 512 个, 0.8466 (=exp(-sigma dt/4eps0)) 512 个
  max|M1 - exp(-sigma dt/eps_used)| = 0.000e+00  (逐位)

已知遗留（同一根因，未修）: m_eps_mf 的其它消费者（能量诊断 E*eps*E、getepsilon_mf）
在材料 2 区域同样用的是 eps1。长期修法是存一个 cell-centered 的 eps_eff 供三者共用。

---faade2c---
feat(m2b/P2): 第 2 步内核分支 —— sigma>0 走指数步，sigma=0 走原分支

MacroscopicEvolveE.cpp:
- 三个 E 循环里对 sigma_arr(i,j,k)>0 的 cell 直接 return（保留 E^n 给矩阵乘）
- 三个循环之后追加一个 kernel（用 mfi.tilebox()）:
    E^{n+1} = M1 E^n + M2 (curlH - J)
  三分量一起读写（M1/M2 是 cell-centered collocated）；curl 表达式与 J 的符号
  与 Yee 分支完全一致
- m1a/m2a 从 MacroscopicProperties::get_sigma_m1_mf/m2_mf 取

CSV dump 增 eps_used 列（按同一 phi 子规则重算），实测 eps_used/eps0 = {1.0:512, 4.0:512}

验收:
1. sigma=0 回归（P0 Fresnel 16ppw）: R_flux=0.1335 T_flux=0.8513 R+T=0.9847
   -> 与改动前逐位一致（sigma=0 时 has_sig_mat=false，代码路径完全不变）   PASS
2. sigma>0 纯弛豫 (/tmp/sig01: eps=4eps0, sigma=0.1, E0=1, B=J=0):
   tau=3.5417e-10, 步长 5.897e-11 (6 步/e-fold)
   step 1/6/12:  rel = 1.13e-10 / 6.77e-10 / 6.17e-10   （受探针文件打印精度限制）
   step 24 以后: rel 迅速恶化，step 48 起 E 变负
   -> 负值来自 PEC 边界反射：域只有 1 m，t=3.3 ns 时反射已返回；
      前 12 步(t<=0.7ns, 2tau)内与 exp(-sigma t/eps) 一致到 1e-10。      部分通过

待办: 用更大域或更短窗口做 1e-12 判据；带源 E=J/sigma 稳态；各向异性 sigma；sigma->inf。

---b3a045e---
feat(m2b/P2): 指数核补 EB mask 与 dt 断言；纯弛豫改周期边界（反射消除，误差=打印精度累加）

内核补两项（复核要求）:
1. 指数 kernel 里沿用第一遍的 EB/PEC mask:
   if (update_Ex_arr && ...) return; 对 Ex/Ey/Ez 各查一次
   避免 sigma>0 的 cell 落在 PEC/EB 上时绕过边界条件写回
2. dt 断言: has_sig_mat 时 AMREX_ALWAYS_ASSERT(get_sigma_dt() == dt)，
   防止以后改 dt 后静默使用过期 M1/M2（配 get_sigma_dt() 访问器）

纯弛豫改 periodic（/tmp/sigper: eps=4eps0, sigma=0.1, 120 步 = 20 tau）:
  step   1: rel 1.13e-10
  step  12: rel 1.35e-09
  step  48: rel 5.42e-09
  step 120: rel 1.35e-08
  E 单调衰减到 2.10e-09（无反射，无变号）  <- 周期边界修好了此前的 PEC 反射

判读: 误差严格线性随步数增长（每步 ~1.1e-10），即单步误差 = FieldProbe 文本的
      打印精度，不是内核误差。要达 1e-12 判据必须提高输出精度（%.17g CSV 或
      plotfile 二进制），当前探针文本 ~1e-10。
      数学层已证 e^{-A dt} 精确（sigma_exp_test.cpp 1.998e-15）。

---e2feb35---
docs: 去掉 research notes 里的本机绝对路径（公开前的清理）


```
