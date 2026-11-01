#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
import numpy as np, sys, os

# columns: 0 step, 1 time, 2 x, 3 y, 4 z, 5 Ex, 6 Ey, 7 Ez, 8 Bx, 9 By, 10 Bz, 11 S
for p in sys.argv[1:]:
    if not os.path.exists(p):
        print(p, "MISSING"); continue
    data = np.loadtxt(p, skiprows=1)
    print(f"{p}: n={len(data)} max|Ex|={np.max(np.abs(data[:,5])):.6e} "
          f"max|Ey|={np.max(np.abs(data[:,6])):.6e} max|Ez|={np.max(np.abs(data[:,7])):.6e}")
