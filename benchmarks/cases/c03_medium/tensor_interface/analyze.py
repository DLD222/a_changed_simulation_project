#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C03-E3: R/T of a tensor-eps interface (normal incidence, s-pol), vacuum-reference method.

  run A: tensor material        diags_A/
  run B: vacuum reference       diags_B/   (no interface, eps = eps0 everywhere)

  R = |E_tot - E_ref|^2 / |E_ref|^2          at the reflection probe
  T = (n2/n1) |E_trans / E_ref|^2            at the transmission probe
  theory: n2^2 = eps2_yy/eps0,  R = ((n1-n2)/(n1+n2))^2,  T = 1 - R

Usage: python3 analyze.py diags_A diags_B [--eps2-yy 5.0]
"""
import argparse
import os
import numpy as np


def load(path):
    d = np.loadtxt(path, skiprows=1)
    return d[:, 1], d[:, 6]          # time, Ey


def amp_at_carrier(t, y, mask):
    yy = y[mask] - y[mask].mean()
    f = np.fft.rfftfreq(len(yy), t[1] - t[0])
    w = 2 * np.pi * f[np.argmax(np.abs(np.fft.rfft(yy)))]
    return np.trapezoid(y[mask] * np.exp(-1j * w * t[mask]), t[mask])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diags_a")
    ap.add_argument("diags_b")
    ap.add_argument("--eps2-yy", type=float, default=5.0)
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--t1", type=float, default=None)
    a = ap.parse_args()

    ta, ya = load(os.path.join(a.diags_a, "qrefl.txt"))
    tb, yb = load(os.path.join(a.diags_b, "qrefl.txt"))
    ta2, ya2 = load(os.path.join(a.diags_a, "qtrans.txt"))
    tb2, yb2 = load(os.path.join(a.diags_b, "qtrans.txt"))
    tmax = min(ta.max(), tb.max()) if a.t1 is None else a.t1
    m = (ta >= a.t0) & (ta <= tmax)
    m2 = (tb >= a.t0) & (tb <= tmax)
    mm = (ta2 >= a.t0) & (ta2 <= tmax)
    mm2 = (tb2 >= a.t0) & (tb2 <= tmax)

    E_ref = amp_at_carrier(tb, yb, m2)
    E_tot = amp_at_carrier(ta, ya, m)
    E_tr = amp_at_carrier(ta2, ya2, mm)
    E_tr_ref = amp_at_carrier(tb2, yb2, mm2)

    n1 = 1.0
    n2 = np.sqrt(a.eps2_yy)
    R = abs(E_tot - E_ref) ** 2 / abs(E_ref) ** 2
    T = (n2 / n1) * abs(E_tr / E_tr_ref) ** 2
    R_th = ((n1 - n2) / (n1 + n2)) ** 2
    T_th = 1 - R_th
    print(f"n2 = {n2:.6f}   R = {R:.6f} (th {R_th:.6f}, {100*(R/R_th-1):+.2f}%)"
          f"   T = {T:.6f} (th {T_th:.6f}, {100*(T/T_th-1):+.2f}%)"
          f"   R+T = {R+T:.6f}")


if __name__ == "__main__":
    main()
