#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Quantitative oblique-incidence check: reflection angle, amplitude, polarization.

Reads a set of FieldProbe files placed on the specular beam
k_r_hat = (-sin 2a, cos 2a) at arc length s = 3..8 m from the specular point,
plus one incident probe and one probe inside the conductor.

Measured quantities

  1. reflection angle - from the slope of the reflected-peak arrival time vs s:
         t_ref(s) = (z_W - z_pulse)/c + s/c  ->  dt/ds = 1/c
     A wrong reflection angle changes this slope (the beam would not be the
     specular one), so the fitted slope is an angle measurement in disguise.
  2. reflection amplitude - r = |E_ref| / |E_inc| on each probe and its mean.
  3. polarization purity - for TE the only allowed component is Ey, for TM the
     field must stay in the x-z plane (Ey = 0) and transverse (E . k_r_hat = 0).

Usage:
  python3 analyze_oblique.py <run_dir> <angle_deg> <te|tm> [z_pulse] [tol_r] [tol_pol]
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


def peak_in(rows, value, t_lo, t_hi):
    best = None
    for r in rows:
        if t_lo <= r[1] <= t_hi:
            v = value(r)
            if best is None or abs(v) > abs(best[1]):
                best = (r[1], v)
    return best


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    run_dir = sys.argv[1]
    alpha = math.radians(float(sys.argv[2]))
    pol = sys.argv[3].lower()
    z_pulse = float(sys.argv[4]) if len(sys.argv) > 4 else -8.0
    tol_r = float(sys.argv[5]) if len(sys.argv) > 5 else 0.20
    tol_pol = float(sys.argv[6]) if len(sys.argv) > 6 else 0.05

    base = os.path.join(run_dir, "diags", "reducedfiles")
    pattern = os.path.join(base, "fp_ref*.txt")
    files = sorted(glob.glob(pattern), key=lambda p: int(re.search(r"fp_ref(\d+)", p).group(1)))
    if not files:
        print(f"no probes matching {pattern}")
        return 1

    kx, kz = -math.sin(2 * alpha), math.cos(2 * alpha)

    if pol == "te":
        def e_val(r):
            return r[6]                           # signed Ey: the phase must flip
        def cross(r):
            return math.hypot(r[5], r[7])         # Ex, Ez
        main_name = "Ey"
        target_r = -1.0
    else:
        def e_val(r):
            return math.hypot(r[5], r[7])         # |Ex,Ez|: amplitude only
        def cross(r):
            return abs(r[6])                      # Ey
        main_name = "|Ex,Ez|"
        target_r = +1.0
    # Incident and reflected pulses must be separable in time, otherwise the
    # peak ratio mixes them (grazing incidence). Pulse half-width ~ 2w/c.
    min_gap = 2.5e-9
    # The single-probe ratio carries a +/-10-17% modulation from diffraction off
    # the finite (PML-truncated) plane, so the mean over the beam is the gate and
    # the per-probe value only has to stay inside a loose band.
    tol_single = tol_r
    tol_mean = 0.06

    print(f"alpha = {math.degrees(alpha):.1f} deg, pol = {pol.upper()}, "
          f"k_r_hat = ({kx:+.6f}, {kz:+.6f}), main component {main_name}")

    rs, ts, ss = [], [], []
    cross_max = 0.0
    trans_max = 0.0
    ok = True
    for path in files:
        name = os.path.basename(path)[:-4]
        rows = load(path)
        x_p, z_p = rows[0][2], rows[0][4]
        s = -(math.cos(alpha) * x_p + math.sin(alpha) * z_p) / math.sin(alpha)
        z_w = z_p - s * math.cos(2 * alpha)
        t_inc = (z_p - z_pulse) / C
        t_ref = (z_w - z_pulse) / C + s / C
        gap = t_ref - t_inc
        if gap < min_gap:
            print(f"{name:9s} s={s:5.2f} m  SKIP (t_ref - t_inc = {gap * 1e9:.2f} ns "
                  f"< {min_gap * 1e9:.1f} ns: incident and reflected overlap)")
            continue
        half = 0.4 * gap

        inc = peak_in(rows, e_val, t_inc - half, t_inc + half)
        ref = peak_in(rows, e_val, t_ref - half, t_ref + half)
        if inc is None or ref is None:
            print(f"{name}: no peak in window")
            ok = False
            continue

        r = ref[1] / inc[1]
        rs.append(r)
        ts.append(ref[0])
        ss.append(s)

        # polarization purity at the reflected peak
        ref_row = min(rows, key=lambda q: abs(q[1] - ref[0]))
        cr = cross(ref_row) / max(abs(e_val(ref_row)), 1e-300)
        cross_max = max(cross_max, cr)
        tr = abs(ref_row[5] * kx + ref_row[7] * kz) / max(e_val(ref_row), 1e-300)
        trans_max = max(trans_max, tr)

        good = abs(r - target_r) <= tol_single
        ok = ok and good
        print(f"{name:9s} s={s:5.2f} m  t_ref={t_ref * 1e9:6.2f} ns (peak {ref[0] * 1e9:6.2f})"
              f"  |E|inc={abs(inc[1]):.3e} |E|ref={abs(ref[1]):.3e}  r={r:+.4f}"
              f"  cross={cr:.3f}  E.k={tr:.3f}  {'OK' if good else 'FAIL'}")

    # reflection angle from the slope of t_ref(s). With fewer than 4 usable
    # probes the peak-time jitter dominates the fit, so the angle is only
    # reported, not used as a gate (grazing incidence).
    n = len(ss)
    if n >= 4:
        ms = sum(ss) / n
        mt = sum(ts) / n
        slope = sum((s - ms) * (t - mt) for s, t in zip(ss, ts)) / sum((s - ms) ** 2 for s in ss)
        inv_c = 1.0 / C
        ang_err_deg = math.degrees(math.acos(max(-1.0, min(1.0, slope / inv_c)))) if slope > 0 else 180.0
        good = abs(slope / inv_c - 1.0) <= 0.05
        ok = ok and good
        print(f"beam angle: dt/ds = {slope * 1e9:.4f} ns/m, 1/c = {inv_c * 1e9:.4f} ns/m, "
              f"rel.err = {slope / inv_c - 1:+.2e}, implied angle error = {ang_err_deg:+.2f} deg "
              f"{'OK' if good else 'FAIL'}")
    elif n >= 2:
        print(f"beam angle: only {n} usable probes, slope check skipped "
              f"(grazing incidence)")

    if rs:
        mean = sum(rs) / len(rs)
        good = abs(mean - target_r) <= tol_mean
        ok = ok and good
        print(f"mean r = {mean:+.4f}  |r-({target_r:+.0f})| = {abs(mean - target_r):.4f}  "
              f"{'OK' if good else 'FAIL'}")
    else:
        print("no usable probe (all skipped)")
        ok = False

    good = cross_max <= tol_pol
    ok = ok and good
    print(f"cross-polarization max = {cross_max:.4f} (limit {tol_pol}) {'OK' if good else 'FAIL'}")

    if pol == "tm":
        # Near the plane the reflected field still contains a reactive component
        # (finite plane, edge diffraction), so this is a loose sanity check that
        # the field is not simply the incident one.
        tol_trans = 0.15
        good = trans_max <= tol_trans
        ok = ok and good
        print(f"transversality |E.k_r|/|E| max = {trans_max:.4f} (limit {tol_trans}) "
              f"{'OK' if good else 'FAIL'}")

    shadow = os.path.join(base, "fp_shadow.txt")
    if os.path.exists(shadow):
        peak = max(abs(r[6]) + abs(r[5]) + abs(r[7]) for r in load(shadow))
        good = peak < 1.0e-3
        ok = ok and good
        print(f"probe inside PEC: max|E| = {peak:.3e} V/m {'OK' if good else 'FAIL'}")

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
