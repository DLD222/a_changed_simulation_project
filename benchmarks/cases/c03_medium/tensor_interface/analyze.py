#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C03-E3: R/T of a tensor-eps interface at normal incidence (vacuum-reference method).

Two changes with respect to the first version (which gave an unstable R):

  1. the carrier is **fixed** to omega = 2*pi*c/lambda (it is an input, not an FFT peak);
  2. each pulse is integrated over its **own** window (centre = emission time + travel
     time, half-width = 5 tau), so the incident and the reflected pulses are never
     truncated by a common window.

  E_inc   = DFT(run B, qrefl, incident window)
  E_refl  = DFT(run A - run B, qrefl, reflected window)
  E_tr    = DFT(run A, qtrans, transmitted window)
  E_trref = DFT(run B, qtrans, transmitted window)
  R = |E_refl|^2 / |E_inc|^2
  T = (n2/n1) |E_tr / E_trref|^2

`--selftest` builds synthetic Gaussian pulses with known R/T and checks that this
analysis recovers them; run it before trusting the WarpX numbers.

Usage:
  python3 analyze.py --selftest
  python3 analyze.py diags_A diags_B --lam 1.0 --t0 16e-9 --tau 2e-9 \
      --z-src -3 --z-refl -2 --z-trans 1 --n2 2.2360679775 --eps2-yy 5.0
"""
import argparse
import os
import numpy as np

C = 299792458.0
EPS0 = 8.8541878128e-12


def load(path):
    d = np.loadtxt(path, skiprows=1)
    return d[:, 1], d[:, 6]          # time, Ey


def dft(t, y, w, t_lo, t_hi):
    m = (t >= t_lo) & (t <= t_hi)
    return np.trapezoid(y[m] * np.exp(-1j * w * t[m]), t[m])


def pulsetimes(t0, tau, z_src, z_refl, z_trans, n2):
    """Pulse centres at the probes (emission time + travel time)."""
    d_inc = (z_refl - z_src) / C
    d_refl = (-z_src + -z_refl) / C                 # source -> interface -> probe
    d_tr_a = (-z_src) / C + z_trans / (C / n2)      # through the interface
    d_tr_b = (z_trans - z_src) / C                  # vacuum reference
    return (t0 + d_inc, t0 + d_refl, t0 + d_tr_a, t0 + d_tr_b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diags_a", nargs="?")
    ap.add_argument("diags_b", nargs="?")
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--t0", type=float, default=16.0e-9)
    ap.add_argument("--tau", type=float, default=2.0e-9)
    ap.add_argument("--z-src", type=float, default=-3.0)
    ap.add_argument("--z-refl", type=float, default=-2.0)
    ap.add_argument("--z-trans", type=float, default=1.0)
    ap.add_argument("--n2", type=float, default=np.sqrt(5.0))
    ap.add_argument("--eps2-yy", type=float, default=5.0)
    ap.add_argument("--nw", type=float, default=5.0, help="window half width in tau")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    w = 2 * np.pi * C / a.lam
    t_inc, t_refl, t_tra, t_trb = pulsetimes(a.t0, a.tau, a.z_src, a.z_refl,
                                             a.z_trans, a.n2)
    hw = a.nw * a.tau

    if a.selftest:
        dt = 3.0e-11
        t = np.arange(0.0, 60.0e-9, dt)
        R_th = ((1 - a.n2) / (1 + a.n2)) ** 2
        T_th = 1 - R_th
        r = np.sqrt(R_th)
        tr = np.sqrt(T_th / a.n2)
        def pulse(tc, amp=1.0, ph=0.0):
            return amp * np.exp(-((t - tc) / a.tau) ** 2) * np.cos(w * (t - tc) + ph)
        qrA = pulse(t_inc) + pulse(t_refl, amp=r, ph=np.pi)
        qrB = pulse(t_inc)
        qtA = pulse(t_tra, amp=tr)
        qtB = pulse(t_trb)
        E_inc = dft(t, qrB, w, t_inc - hw, t_inc + hw)
        E_refl = dft(t, qrA - qrB, w, t_refl - hw, t_refl + hw)
        E_tr = dft(t, qtA, w, t_tra - hw, t_tra + hw)
        E_trref = dft(t, qtB, w, t_trb - hw, t_trb + hw)
        R = abs(E_refl) ** 2 / abs(E_inc) ** 2
        T = a.n2 * abs(E_tr / E_trref) ** 2
        print(f"selftest: R = {R:.6f} (th {R_th:.6f})  T = {T:.6f} (th {T_th:.6f})  "
              f"R+T = {R + T:.6f}")
        ok = abs(R / R_th - 1) < 1e-6 and abs(T / T_th - 1) < 1e-6
        print("SELFTEST PASS" if ok else "SELFTEST FAIL")
        return 0 if ok else 1

    ta, ya = load(os.path.join(a.diags_a, "qrefl.txt"))
    tb, yb = load(os.path.join(a.diags_b, "qrefl.txt"))
    ta2, ya2 = load(os.path.join(a.diags_a, "qtrans.txt"))
    tb2, yb2 = load(os.path.join(a.diags_b, "qtrans.txt"))

    E_inc = dft(tb, yb, w, t_inc - hw, t_inc + hw)
    E_refl = dft(ta, ya, w, t_refl - hw, t_refl + hw) - \
             dft(tb, yb, w, t_refl - hw, t_refl + hw)
    E_tr = dft(ta2, ya2, w, t_tra - hw, t_tra + hw)
    E_trref = dft(tb2, yb2, w, t_trb - hw, t_trb + hw)
    R = abs(E_refl) ** 2 / abs(E_inc) ** 2
    T = a.n2 * abs(E_tr / E_trref) ** 2
    R_th = ((1 - a.n2) / (1 + a.n2)) ** 2
    T_th = 1 - R_th
    print(f"n2 = {a.n2:.6f}")
    print(f"R = {R:.6f}  th {R_th:.6f}  dev {100*(R/R_th-1):+.2f}%")
    print(f"T = {T:.6f}  th {T_th:.6f}  dev {100*(T/T_th-1):+.2f}%")
    print(f"R+T = {R+T:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
