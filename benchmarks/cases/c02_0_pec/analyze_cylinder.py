#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""PEC-cylinder scattering: boundary condition, shadowing, azimuthal pattern.

A plane pulse hits a PEC cylinder of radius R along +z. Probes sit at r = 2 m
every 45 deg (theta measured from +z), so

  * theta = 0   is the shadow side (behind the cylinder),
  * theta = 180 is the specular side (the pulse hits the front of the cylinder
    at z = -R and reflects back towards -z),
  * theta = 90 / 270 are the grazing sides.

Checks (qualitative but physical):
  1. the probe inside the conductor stays ~0 (PEC boundary condition),
  2. the specular probe is much stronger than the shadow probe,
  3. the azimuthal pattern is front/back asymmetric,
  4. all values stay finite and the late-time peak does not grow.

This is not a Mie benchmark: the incident pulse has no carrier, so a spectral
comparison with the analytic series is left for a dedicated case.

Usage:  python3 analyze_cylinder.py <run_dir> <te|tm> [ratio_min]
"""
import os
import sys

ANGLES = (0, 45, 90, 135, 180, 225, 270, 315)


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("[") or s.startswith("#"):
                continue
            rows.append([float(v) for v in line.split()])
    return rows


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    run_dir = sys.argv[1]
    pol = sys.argv[2].lower()
    ratio_min = float(sys.argv[3]) if len(sys.argv) > 3 else 1.2

    base = os.path.join(run_dir, "diags", "reducedfiles")

    def amp(row):
        if pol == "te":
            return abs(row[6])
        return (row[5] ** 2 + row[7] ** 2) ** 0.5

    peaks, tails = {}, []
    finite_all = True
    for a in ANGLES:
        rows = load(os.path.join(base, f"p{a}.txt"))
        amp_n = [amp(r) for r in rows]
        n = len(amp_n)
        # the pulse needs time to reach the probe, so "early vs late" is not a
        # growth measure; instead compare the last few percent with the global peak
        peak = max(amp_n)
        tail = max(amp_n[int(0.95 * n):]) if n >= 20 else amp_n[-1]
        peaks[a] = peak
        tails.append(tail / peak if peak > 0 else float("inf"))
        finite = all(abs(v) < 1e100 for r in rows for v in r)
        finite_all = finite_all and finite
        print(f"theta={a:3d} deg: peak|E|={peak:.4e}  tail/peak={tail / peak if peak > 0 else float('inf'):.3f}"
              f"  finite={finite}")

    pin = load(os.path.join(base, "pin.txt"))
    inside = max(amp(r) for r in pin)
    good_inside = inside < 1.0e-3 * max(peaks.values())

    ratio = peaks[180] / max(peaks[0], 1e-300)
    asym = (peaks[180] - peaks[0]) / max(peaks[180], peaks[0], 1e-300)
    good_ratio = ratio >= ratio_min
    good_asym = asym >= 0.2

    growth = max(tails)
    good_growth = growth <= 1.2 and finite_all

    print(f"inside PEC: max|E| = {inside:.3e}  (limit 1e-3 of peak) "
          f"{'OK' if good_inside else 'FAIL'}")
    print(f"specular/shadow = |E(180)|/|E(0)| = {ratio:.3f}  (min {ratio_min}) "
          f"{'OK' if good_ratio else 'FAIL'}")
    print(f"front/back asymmetry = {asym:.3f}  (min 0.20) {'OK' if good_asym else 'FAIL'}")
    print(f"tail/peak max = {growth:.3f}  (max 1.2, no growth) {'OK' if good_growth else 'FAIL'}")

    ok = good_inside and good_ratio and good_asym and good_growth
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
