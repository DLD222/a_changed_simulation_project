# 介质相速度判据修正、根因定位（参数实验）与修改意见

> 日期：2026-10-31 · 对象：`MEDIUM_STATUS_REPORT.md` 所述"均匀介质三项 FAIL"问题
> 方法：**跑算例 + 逐参数做实验**（不修改原报告），用对照实验把问题从"材料系数/内核"逼到"包装函数"；
> 再用**宽脉冲 + 包络峰值**的测速判据验证修复。
> 结论：根因是 `WarpX::MacroscopicEvolveE` 漏了 `ECTRhofield` 重算；修复后
> **相速度误差 0.1–0.3%、介质加权能量守恒 1e-4**。本文最后给出**修改意见（建议性）**。

---

## 1. 旧判据为什么不可靠

| 问题 | 说明 |
|---|---|
| 脉冲只有 3.6 格 | 初始高斯宽 0.15 m、dx=4/96≈0.042 m；窄脉冲数值色散极强，"首达时间"由最短波长决定，必然偏快 |
| 用"首个 5% 穿越" | 色散展宽的脉冲前沿到达时间 ≠ 包络峰值时间 |
| 真空对照也报 +203% | ECT 真空（物理 v=c）被测成约 2c——判据本身失效 |
| 初始 B 未按介质匹配 | 原测试给 `B=−E/c`，在 εr≠1 里是 ± 行波混合，驻波节点让探测点接近 0 |

## 2. 修正判据（宽脉冲 + 包络峰值到达时间）

- 高斯包络 `w = 0.5 m`（128² 时 ≈12 格），载波 `λ0 = 1 m`（24 ppw）；
- `Ey = A(z)cos(k0(z−z0))`、`Bx = −Ey/v`（按介质匹配，单向 +z）；
- 三探针 z = −0.5 / 0 / 0.5，取**前半段**（反射返回前）`|Ey|` 峰值时刻，`v = Δz/Δt`；
- 能量用**介质加权** `εr·E_lev0 + B_lev0`（`FieldEnergy` 原始输出是 ε0 加权，介质里不是守恒量）；
- 脚本：`benchmarks/cases/c02_0_pec/verify_medium_speed.py`。

## 3. 结果（修复后）

| εr | v_theory | v_peak（p1→p3） | 相对误差 | 能量（介质加权）final/init |
|---|---|---|---|---|
| 1.00 | 1.0000c | 1.0012c | **+0.12%** | 0.99991 |
| 1.50 | 0.8165c | 0.8139c | **−0.31%** | 0.99993 |
| 2.25 | 0.6667c | 0.6655c | **−0.17%** | 0.99997 |

对比：旧判据在 ECT 真空上就报 `v≈2c`（+203%），修正后三种 εr 均在 **0.3%** 内。

## 4. 根因定位：参数实验矩阵

按"每步只改一个参数/一行代码"做对照（同一算例、同一网格、同一二进制流程）：

| # | 实验（改动） | 能量比值 | 判定 |
|---|---|---|---|
| 1 | 原样（εr=2.25，ECT+材料） | ×1.663e5 | 复现报告 |
| 2 | **εr=1.0**（材料路径不变） | ×**8.401e5** | 与 ε **无关** |
| 3 | 临时令 `mur_inv=1`（去掉 μ 因子） | ×8.401e5 不变 | 与 μ **无关** |
| 4 | 临时把材料分支换成纯 `E += dE` | ×8.401e5 不变 | 与材料系数 **无关** |
| 5 | 临时强制 `has_macro=false`（内核=真空） | ×8.401e5 不变 | **与 ECT 内核无关** |
| 6 | 内核调用计数 | 每步 **1 次** | "E 被更新两次"（报告假设 A）**排除** |
| 7 | 对比 `WarpX::EvolveE` 与 `WarpX::MacroscopicEvolveE` 源码 | — | 真空路径末尾有 `EvolveECTRho` 重算，**宏观路径没有** |
| 8 | 补上 `EvolveECTRho` 后重跑 | εr=1 **0.9999**；εr=2.25 介质加权 **1.0000（0.04%）** | **根因确认** |

补充观察（与根因一致）：
- E 能量在第 1 步与真空**完全一致**，但 **B 能量第 1 步就差 0.6%**——B 用的正是陈旧的 `ECTRhofield`；
- 增长不是 CFL 型（cfl=0.4/0.2/0.1 的每步增长率 3.4%/2.7%/2.0%，不随 dt 线性下降）；
- 能量由 E 主导（末态 E_lev0=5.5e-6、B_lev0=4.1e-9）。

## 5. 根因

`WarpX::MacroscopicEvolveE`（`Source/FieldSolver/WarpXPushFieldsEM.cpp`，`(lev, PatchType, dt, time)` 重载）
在 E 更新后只调用 `ApplyEfieldBoundary`，**漏掉了**真空路径末尾的 `ECTRhofield` 重算：

```cpp
// WarpX::EvolveE 里有，MacroscopicEvolveE 里没有：
if (electromagnetic_solver_id == ECT) { ... EvolveECTRho(Efield, edge_lengths, face_areas, ECTRhofield, lev); }
```

于是 ECT 的 B 更新（`EvolveBCartesianECT` 用 `ECTRhofield`）读的是 **E 更新之前**的 Rho，
每步注入能量并形成正反馈 → 指数增长。这解释了：εr=1 也爆、μ/材料系数/内核都无关、
真空路径正常、Yee+介质正常（Yee 路径不用 `ECTRhofield`）。

---

## 6. 修改意见（建议性，按优先级）

### A. 代码修复（根因，已验证）

在 `WarpX::MacroscopicEvolveE` 的 `ApplyEfieldBoundary(lev, patch_type, new_time);` 之后，
补上与 `WarpX::EvolveE` 相同的 `EvolveECTRho` 调用（fine/coarse 各一份，`#ifdef AMREX_USE_EB`
+ `solver_id == ECT` 守卫）：

```cpp
#ifdef AMREX_USE_EB
    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT) {
        if (patch_type == PatchType::fine) {
            m_fdtd_solver_fp[lev]->EvolveECTRho(
                m_fields.get_alldirs(FieldType::Efield_fp, lev),
                m_fields.get_alldirs(FieldType::edge_lengths, lev),
                m_fields.get_alldirs(FieldType::face_areas, lev),
                m_fields.get_alldirs(FieldType::ECTRhofield, lev), lev );
        } else { /* 同款，用 Efield_cp / solver_cp */ }
    }
#endif
```

建议同时做**防再犯**：把"E 更新后重算 ECTRhofield"抽成一个共享函数（或 `WarpX` 私有方法），
`EvolveE` 与 `MacroscopicEvolveE` 都调用它；根因正是这段代码被复制后被漏掉。

### B. 测试/诊断修复（让这类 bug 能被抓住）

1. `verify_medium_wave.py` 的能量判据改**介质加权**：`εr·E_lev0 + B_lev0`（或同时输出 E、B）；
2. 速度判据改**宽脉冲包络峰值**（或互相关），去掉"5% 首达"；脉冲宽 ≥10 格、`B=−E/v`；
3. 增加 εr=1.5 一档；**保留 εr=1 作为 sanity check**（修复前它也会爆，正好抓本 bug）；
4. 回归里加"ECT+介质 + `ect_conformal_E>0`"的组合（当前测试没有覆盖 ECT×介质）。

### C. 健壮性建议

1. 加一个**调试断言**：在 ECT+macroscopic 的 E 更新末尾，对比"重算前/后"的 `ECTRhofield`
   （或检查其与 `EvolveECTRho` 结果的差异），不一致即告警；
2. 在 `MacroscopicEvolveE` 的入口/出口加一次 `WarpX::verbose` 级打印（路径选择），便于今后定位；
3. 参考 `MEDIUM_CODE_AUDIT.md` 的 P2/P3/P4/P6/P8/P9 排期：P6（界面回归）与本 bug 同属"测试缺口"，
   建议一起补。

### D. 与 `MEDIUM_STATUS_REPORT.md` 的关系

- 该报告的三项 FAIL 数字（1.663e5 / −27.7% / B/A=0）**复现无误**；
- 其"两条待验证假设"（双更新、非线性反馈）经参数实验**均排除**；根因是本文 §5；
- 其"介质不可用、PEC 路径不受影响"的结论仍成立；
- 建议把本文件作为报告附录/勘误引用（是否并入由项目决定，本文不修改原报告）。

---

## 7. 验证清单（目标）

| 判据 | 目标 | 实测（修复后） |
|---|---|---|
| 相速度 | \|v/v_th−1\| ≤ 1% | 0.12% / 0.31% / 0.17%（εr=1/1.5/2.25） |
| 能量（介质加权） | 漂移 ≤ 1e-4 | ≤ 1e-4 |
| 振幅比 B/A | 0.5–2.0 | 0.977（旧判据口径） |
| 界面（Fresnel） | \|R\| ≤1% | 待做（M2b V1） |

## 8. 复现

```bash
# 修复后的测速（宽脉冲 + 包络峰值）
python3 benchmarks/cases/c02_0_pec/verify_medium_speed.py /tmp/spd_225 --epsr 2.25 --make --execute
# 参数实验矩阵：#2 εr=1、#3 μ 关、#4 E+=dE、#5 has_macro=false、#6 调用计数
# 见本文 §4（均为临时改动，验证后已还原）
```

*本文为建议性记录：代码修复已在工作区（`WarpX::MacroscopicEvolveE`，未提交）；是否按 A–D 推进由项目决定。*
</content>
