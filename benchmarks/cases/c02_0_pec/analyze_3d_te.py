#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""3D TE pulse check: propagation speed, amplitude, Bx/Ey ratio.

Reads two point probes on the z axis (p1 at z1, p2 at z2 > z1) of a planar TE
pulse initialized at z = z_pulse with Bx = -Ey/c. For a wave travelling at c:

  * t(p2) - t(p1) = (z2 - z1)/c,
  * the peak |Ey| is the same at both probes (planar, no geometric spreading),
  * Bx/Ey = -1/c.

Usage:  python3 analyze_3d_te.py <run_dir> [z_pulse] [tol_time_ns] [tol_amp]
"""
import os
import sys

C = 299792458.0


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("[") or s.startswith("#"):
                continue
            rows.append([float(v) for v in line.split()])
    return rows


def peak_in(rows, column, t_lo, t_hi):
    best = None
    for r in rows:
        if t_lo <= r[1] <= t_hi and (best is None or abs(r[column]) > abs(best[1])):
            best = (r[1], r[column])
    return best


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    run_dir = sys.argv[1]
    z_pulse = float(sys.argv[2]) if len(sys.argv) > 2 else -2.0
    tol_time_ns = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    tol_amp = float(sys.argv[4]) if len(sys.argv) > 4 else 0.10

    base = os.path.join(run_dir, "diags", "reducedfiles")
    data = {}
    for name in ("p1", "p2"):
        rows = load(os.path.join(base, name + ".txt"))
        z = rows[0][4]
        t_exact = (z - z_pulse) / C
        half = 0.5 * abs(t_exact - (rows[0][4] - z_pulse) / C) + 6.0e-9
        window = (t_exact - half, t_exact + half)
        ey = peak_in(rows, 6, *window)
        bx = peak_in(rows, 8, *window)
        if ey is None or bx is None:
            print(f"{name}: no peak found near {t_exact * 1e9:.2f} ns")
            return 1
        data[name] = {"z": z, "t": ey[0], "ey": ey[1], "bx": bx[1], "t_exact": t_exact}
        print(f"{name}: z={z:+.3f} m  t_peak={ey[0] * 1e9:7.3f} ns "
              f"(analytic {(z - z_pulse) / C * 1e9:7.3f})  Ey={ey[1]:+.4e}  Bx={bx[1]:+.4e}")

    dz = data["p2"]["z"] - data["p1"]["z"]
    dt_exact = dz / C
    dt = data["p2"]["t"] - data["p1"]["t"]
    amp_ratio = abs(data["p2"]["ey"]) / abs(data["p1"]["ey"])
    bx_ey = data["p1"]["bx"] / data["p1"]["ey"]
    bx_ey_exact = -1.0 / C

    print(f"propagation: dt_measured={dt * 1e9:.3f} ns  dt_analytic={dt_exact * 1e9:.3f} ns  "
          f"err={(dt - dt_exact) * 1e9:+.3f} ns")
    print(f"amplitude:   |Ey(p2)|/|Ey(p1)| = {amp_ratio:.4f}")
    print(f"impedance:   Bx/Ey = {bx_ey:+.4e}  (analytic {bx_ey_exact:+.4e}, "
          f"relerr {abs(bx_ey / bx_ey_exact - 1):.2e})")

    # Bx and Ey live half a time step apart on the Yee grid, so the ratio taken
    # at the same sample is off by O(omega*dt/2) - a few percent for this pulse
    # and grid. The tolerance is therefore loose; the ratio is a sanity check on
    # the sign and magnitude, not a precision measurement.
    checks = [
        (f"|dt - dz/c| = {abs(dt - dt_exact) * 1e9:.3f} ns <= {tol_time_ns} ns",
         abs(dt - dt_exact) <= tol_time_ns * 1e-9),
        (f"|amplitude ratio - 1| = {abs(amp_ratio - 1):.4f} <= {tol_amp}",
         abs(amp_ratio - 1.0) <= tol_amp),
        (f"|Bx/Ey + 1/c| relative = {abs(bx_ey / bx_ey_exact - 1):.2e} <= 0.1 "
         f"(half-step staggering)", abs(bx_ey / bx_ey_exact - 1.0) <= 0.1),
    ]
    ok = True
    for label, good in checks:
        ok = ok and good
        print(f"  {'OK  ' if good else 'FAIL'} {label}")

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
