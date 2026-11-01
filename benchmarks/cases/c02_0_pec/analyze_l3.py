#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""L3: reflection coefficient of a slanted PEC plane, from probes ON the beam.

Geometry: eb = cos(alpha)*(x-x0) + sin(alpha)*(z-z0), PEC on the positive side,
incident pulse travelling along +z. Specular reflection gives

    k_r_hat = (-sin(2 alpha), cos(2 alpha))        (= (-0.866, 0.5) for 30 deg)

For a probe at P the specular point W on the plane follows from
n.P - s (n.k_r_hat) = 0 with n = (cos alpha, sin alpha), n.k_r_hat = -sin(alpha),
hence s = -(n.P)/sin(alpha), z_W = z_P - s*cos(2 alpha), and

    t_inc = (z_P - z_pulse)/c ,   t_ref = (z_W - z_pulse)/c + s/c

Both peaks are searched in a window around those analytic times. A PEC reflector
must give r = Ey_ref/Ey_inc = -1: same amplitude, inverted phase.

The single-probe value is modulated by diffraction off the finite plane and by
the absorbing outer boundary (measured: -0.92 .. -1.11 over s = 3..8 m, and
essentially unchanged between 384^2 and 768^2), so the criterion is

  * every probe shows phase inversion (r < 0),
  * the mean over the beam satisfies |r_mean + 1| <= tol_mean (default 0.05),
  * no probe deviates by more than tol_single (default 0.20),
  * the probe inside the conductor stays negligible.

Usage:  python3 analyze_l3.py <run_dir> [alpha_deg] [z_pulse] [tol_mean] [tol_single]
"""
import glob
import math
import os
import re
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
    alpha = math.radians(float(sys.argv[2]) if len(sys.argv) > 2 else 30.0)
    z_pulse = float(sys.argv[3]) if len(sys.argv) > 3 else -8.0
    tol_mean = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    tol_single = float(sys.argv[5]) if len(sys.argv) > 5 else 0.20

    pattern = os.path.join(run_dir, "diags", "reducedfiles", "fp_ref*.txt")
    files = sorted(glob.glob(pattern), key=lambda p: int(re.search(r"fp_ref(\d+)", p).group(1)))
    if not files:
        print(f"no probes matching {pattern}")
        return 1

    print(f"alpha = {math.degrees(alpha):.1f} deg, k_r_hat = "
          f"({-math.sin(2 * alpha):+.6f}, {math.cos(2 * alpha):+.6f})")
    print(f"criterion: all r < 0, |mean(r)+1| <= {tol_mean}, "
          f"|r+1| <= {tol_single} per probe")

    ok = True
    rs = []
    for path in files:
        name = os.path.basename(path)[:-4]
        rows = load(path)
        x_p, z_p = rows[0][2], rows[0][4]
        s = -(math.cos(alpha) * x_p + math.sin(alpha) * z_p) / math.sin(alpha)
        z_w = z_p - s * math.cos(2 * alpha)
        t_inc = (z_p - z_pulse) / C
        t_ref = (z_w - z_pulse) / C + s / C
        half = 0.45 * (t_ref - t_inc)

        inc = peak_in(rows, 6, t_inc - half, t_inc + half)
        ref = peak_in(rows, 6, t_ref - half, t_ref + half)
        if inc is None or ref is None:
            print(f"{name}: no peak in window ({t_inc * 1e9:.2f} / {t_ref * 1e9:.2f} ns)")
            ok = False
            continue
        r = ref[1] / inc[1]
        rs.append(r)
        good = (r < 0.0) and (abs(r + 1.0) <= tol_single)
        ok = ok and good
        print(f"{name:9s} s={s:5.2f} m  t_inc={t_inc * 1e9:6.2f} ns (Ey={inc[1]:+.4e})"
              f"  t_ref={t_ref * 1e9:6.2f} ns (Ey={ref[1]:+.4e})"
              f"  r={r:+.4f}  {'OK' if good else 'FAIL'}")

    shadow = os.path.join(run_dir, "diags", "reducedfiles", "fp_shadow.txt")
    if os.path.exists(shadow):
        peak = max(abs(r[6]) for r in load(shadow))
        good = peak < 1.0e-3
        ok = ok and good
        print(f"fp_shadow (inside PEC): max|Ey| = {peak:.3e} V/m {'OK' if good else 'FAIL'}")

    if rs:
        mean = sum(rs) / len(rs)
        good = abs(mean + 1.0) <= tol_mean
        ok = ok and good
        print(f"mean over beam: r = {mean:+.4f}  |r+1| = {abs(mean + 1.0):.4f}  "
              f"{'OK' if good else 'FAIL'}")

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
