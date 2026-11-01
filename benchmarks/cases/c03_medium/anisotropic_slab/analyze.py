#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C03-E4: bulk tensor phase speed from the envelope-peak transit time.

  v = dz / (t_peak(q2) - t_peak(q1)),   v_theory = c / sqrt(eps_ii/eps0)
  s-pol (Ey): eps_yy = 5   p-pol (Ex): eps_xx = 4

Usage: python3 analyze.py diags/reducedfiles --comp Ey --eps-ii 5.0 --dz 3.0
"""
import argparse
import os
import numpy as np

C = 299792458.0
COLS = {"Ex": 5, "Ey": 6, "Ez": 7}


def load(path, comp):
    d = np.loadtxt(path, skiprows=1)
    return d[:, 1], d[:, COLS[comp]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diag", nargs="?", default="diags/reducedfiles")
    ap.add_argument("--comp", default="Ey")
    ap.add_argument("--eps-ii", type=float, default=5.0)
    ap.add_argument("--dz", type=float, default=3.0)
    ap.add_argument("--t0", type=float, default=0.0)
    a = ap.parse_args()

    t1, y1 = load(os.path.join(a.diag, "q1.txt"), a.comp)
    t2, y2 = load(os.path.join(a.diag, "q2.txt"), a.comp)
    m1 = t1 >= a.t0
    m2 = t2 >= a.t0
    i1 = np.where(m1)[0][np.argmax(np.abs(y1[m1]))]
    i2 = np.where(m2)[0][np.argmax(np.abs(y2[m2]))]
    dt = t2[i2] - t1[i1]
    v = a.dz / dt
    v_th = C / np.sqrt(a.eps_ii)
    print(f"peak times: {t1[i1]:.6e} / {t2[i2]:.6e} s   dt = {dt:.6e} s")
    print(f"v = {v:.6e} m/s   v_th = {v_th:.6e} m/s   dev {100*(v/v_th-1):+.3f}%")


if __name__ == "__main__":
    main()
