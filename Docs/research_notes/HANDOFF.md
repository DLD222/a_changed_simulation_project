# 交接摘要（2026-10-31 收尾）

<!-- SPDX-License-Identifier: BSD-3-Clause -->


> 目的：把当前状态压缩成一页，后续任何会话从这里开始即可，不必重读全部 md。
> 权威细节见文末索引；本文只保留"结论 + 数字 + 待办"。

## 1. 项目与代码

- 目标：WarpX 分支 `~/WarpX_fluid` 做 **n 组分矩流体 + 完整 Maxwell + 真 cut-cell(EB) + 介质 + AMR**；
  当前攻坚 **M2a = PEC 上 E+B 全分量共形 cut-cell**（2D XZ 优先，3D 同样要求）。
- 构建：`build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB`、`build/bin/warpx.3d...EB`（EB=ON、MPI、双精度）。
- 相对上游（`~/WarpX_changed`）改了 **6 个源文件**：`EmbeddedBoundaryInit.cpp`、`WarpXFaceExtensions.cpp`、
  `WarpXFaceInfoBox.H`、`EvolveB.cpp`、`EvolveECTRho.cpp`、`WarpXInitData.cpp`。
- 关键修复：`EvolveECTRho` 从 `InitLevelData` 移到 `AddExternalFields` 末尾（原调用早于初始 E 注入，
  导致 `ρ≡0`、首半步 B 空转、整体滞后半拍）。2D `idim<3` 让 Bz 参与；XZ 的 `Rhox/Rhoz` 放开；
  per-unit-y 几何（`ly`=节点开放度、`face_areas[0]/[2]`=`getAreaFrac×dz/dx`）；TE 面 BCK-only。
- 调试残留已清理（`DBGB`/`DBGRHO` 已删）。

## 2. 判据状态

| 判据 | 内容 | 结果 |
|---|---|---|
| L1 | 几何退化（`Sx={0,dz}`、`Sz={0,dx}`、`ly={0,1}`） | PASS |
| L2a | 规则域（EB 移出域外）ECT vs Yee | **201 步 `max\|d\|=0`**（逐位） |
| L2b | 网格对齐 PEC 的等效壁面 | ECT 往返误差 `+0.371/+0.017/+0.076 ns`（64/128/256），Yee `−1.75/−0.81/−0.34`；PASS |
| L3 | 斜入射（反射束 6 探针 + 斜率测角） | 45°/60° × TE/TM 四组 PASS（平均 `\|r∓1\|=0.026–0.051`，交叉偏振 0，斜率 0.1–2.4%）；**30° 掠入射 FAIL**（斜率 10.5%，TM 横向性 0.446） |
| L4 | PEC 圆柱（边界/阴影/方位） | PASS：柱内 0；镜面/阴影 TE 5.35 / TM 1.42 |
| 长时间 | 45° 斜面 10000 步 | TE 能量漂移 `3.0e-4`；**TM 1.9% 增长**（未发散，观察项） |
| 上游回归 | 2D/3D `rotated_cube` | `6.8e-10` / `≤1e-13` |
| TM 基线 | `inputs_2d_tm_flat` n=384 | `r = −0.999361` |

## 3. 阻塞项（最重要）

**圆柱米氏收敛阶测不出来**：80²/160² 数十步内爆到 `1e306–1e308`。

**原结论（"根因是 E 侧阶梯，正解 = C 侧"）已于 2026-10-31 晚被实验推翻**，见 `C_SIDE_IMPL_LOG.md` §3：

| 配置（80², 400 步） | step 200 的 `max\|E\|` |
|---|---|
| ECT + EB + 共形 E | 1e308 |
| ECT + EB 移出域外 + 共形 E | 1.36e308 |
| **Yee + EB 移出域外 + PML** | **1.67e308** |
| **Yee + EB 移出域外 + PEC** | **1.46e308** |
| 同上但源幅值置零 | 2.88e6（稳定） |

⇒ 不用 EB、不用 ECT、只用 Yee 也爆；源置零就稳定。**根因在 CW 激光源
（`parse_field_function` + `cos(ωt)`），与 EB/ECT/C 侧无关**，C 侧修不了它。
下一步先定位源（换 `sin`、平滑开启、扫描 40²/64²/80²/128²）。

## 4. 待办（按优先级）

1. ~~**C 侧实现**~~ **已完成**（tag `c-side-2d-v1`）：`EvolveECartesianECT` + `dual_A_y`/`dual_L_x`/`dual_L_z`，
   开关 `algo.ect_conformal_E`（0 关 / 1 常数几何 / 2 volfrac 几何）。
   实测：L2a 到机器精度（`rel ≤ 2e-15`，**不是逐位 0**）；L3 45° 斜率 −2.39% → −1.51%；
   米氏收敛阶近场 1.28（掩膜 Yee）vs 1.06（共形 E），RCS 1.76 / 1.77。
   **仍未做**：掩膜改 `A_y<=0`、对偶稳定化 `A_stab`、方案 A/C 几何、3D、AMR。
   修正清单见 `C_SIDE_FIX_PLAN.md`（F1–F7 + 6 点微调）。
2. **米氏收敛阶**：**已解锁**。原"爆炸"根因是 `max_grid_size = 64 < n_cell`（多 grid + CW 源），
   改成单 grid 后 40²/80²/160² 全部稳定。TE 三档已测（见上一条）。待做：4 档以上（加 320）、
   `RCS` 减**数值入射**（`C_SIDE_FIX_PLAN.md` F6）、真空残差作为误差地板并列报告。
3. **AMR 多层几何**：几何目前只在 `maxLevel()` 初始化；需 per-level 初始化 + 扩展链 + 粗-细一致化。
4. **2D/3D 逐点一致性**：1 层 y + 周期 + `dy=1` + 同 `dx/dz/dt`，容差待定（cut cell 下不保证 1%）。
5. **TM 1.9% 能量增长**：C 侧完成后复测；若仍在，查 `Rhoy` 与 BCK 的交互。
6. 小项：`EvolveB` 借用索引（`idim=0/2` 仍是 3D-y，当前不可达）、`S` 就地覆盖、`ly` 判据偏保守、
   初始化覆盖仅 external-E 路径。

## 5. 关键源码事实（避免重复踩）

- `eb_implicit_function` **正值为固体**（写反会把外部全变导体）。
- 网格数必须是 `blocking_factor`(8) 的倍数。
- `yee + EB + PML` 会 abort（投影 div(B) 清洁器只接受 periodic/PEC/PMC/neumann）；圆柱必须用 ECT。
- `probe_geometry = Line` 在本 2D 构建下返回全 0，不可用于空间定位（用 Full plotfile）。
- `laser1.profile = parse_field_function` 是**可用的 CW 平面波源**：`field_function(X,Y,t)` 变量为
  `{X,Y,t}`，返回值为场幅值(V/m)；2D ECT + TE 偏振已验证。
- 对齐 PEC 上 AMReX 把"贴面单元"标为 cut（`volfrac=1` 但 `isRegular()=false`），Yee 的 staircase 掩膜
  因此多冻结一层 ⇒ ECT 与 Yee 必然会差一层 `dz`，这是 L2b 的全部来源，不是 bug。

## 6. 文件索引

**文档**（`Docs/research_notes/`）
- `ECT_2D_TE_REPORT.md`：主报告（§2 判据表、§4 改动清单、§14 迭代记录）
- `CASES_2026-10-31.md`：**所有算例的物理情景/结果/问题**（C1–C15、§8 八条坑）
- `M2A_CODE_AUDIT.md`：第三方审计 + §0 处置状态
- `L3_REFLECTION_CRITERION.md`：30°/60°/15° 判据
- `C_SIDE_DESIGN.md`：C 侧实现清单 + AMR/2D-3D 方案（**§3–§5 已被 DESIGN2 取代**）
- `C_SIDE_RESEARCH.md`：C 侧**文献调研**（教科书/论文/开源）+ 可修复性判定 + S0–S10 实施流程
- `C_SIDE_DESIGN2.md`：**C 侧实现权威版**（修正公式、对偶几何取法、S0'–S11'、自检 G-1/G-2/F-1/F-2/E-1）
- `C_SIDE_RESEARCH2.md`：C 侧调研第二版（DESIGN2 逐条核对、PML/半格等 P1、16 条存疑）
- `C_SIDE_IMPL_LOG.md`：**C 侧实施记录**（S0'–S3' 代码清单、L2a/L3 结果、米氏爆炸归因的对照实验）
- `A_TM_BASELINE.md`、`B_2D_TE_RESEARCH.md`、`ROOTCAUSE_ECT_2D.md`、`M2A_SOURCE_REVIEW_FINDINGS.md`、
  `cutcell-research.md`、`architecture-final.md`、`design-plan-moment-solver.md`

**算例/脚本**（`benchmarks/cases/c02_0_pec/`）
- 判据：`run_l2.sh`（L2a/L2b）、`run_l2b_scan.sh`（分辨率）、`verify_l2a.py`、`verify_l2b.py`
- 斜入射：`make_oblique.py`、`analyze_oblique.py`、`inputs_2d_oblique_{te,tm}_{30,45,60}deg`
- 圆柱：`make_cylinder.py`、`analyze_cylinder.py`、`inputs_2d_cyl_{te,tm}`
- 米氏：`make_mie_cases.py`、`mie_cylinder.py`、`analyze_mie_convergence.py`、`inputs_mie_{ey,by}_{40..640}`
- 稳定性：`inputs_ect_2d_45deg[_tm]`、`check_stability.py`
- 3D：`inputs_3d_te_pulse`、`analyze_3d_te.py`

## 7. 一键复现

```bash
cd <repo>/benchmarks/cases/c02_0_pec
./run_l2.sh && ./run_l2b_scan.sh          # L2a/L2b + 分辨率扫描
python3 make_oblique.py 45 tm && python3 analyze_oblique.py <run> 45 tm
python3 make_cylinder.py te && python3 analyze_cylinder.py <run> te
python3 mie_cylinder.py                    # Mie 级数自检
python3 make_mie_cases.py                  # 生成米氏收敛阶算例（当前会爆，见 §3）
```
