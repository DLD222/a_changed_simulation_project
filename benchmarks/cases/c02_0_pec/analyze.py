#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C2-0 analysis: extract reflection coefficient from FieldProbe point probes.
columns: 0 step, 1 time, 2 x, 3 y, 4 z, 5 Ex, 6 Ey, 7 Ez, 8 Bx, 9 By, 10 Bz, 11 S"""
import numpy as np

def load(path):
    data = np.loadtxt(path, skiprows=1)
    return data[:,1], data[:,6]   # time, Ey (s-polarization)

t, Ey = load("diags/reducedfiles/fp_inc.txt")
i_peak = np.argmax(np.abs(Ey))
sign = np.sign(Ey[i_peak])
# reflected pulse: opposite-sign extremum strictly after the incident peak
late = (t > t[i_peak] + 0.2 * (t.max() - t.min()))
opp = late & (np.sign(Ey) == -sign)
print(f"incident peak: t={t[i_peak]:.4e} s  Ey={Ey[i_peak]:.6e} V/m")
if opp.any():
    j = np.where(opp)[0][np.argmax(np.abs(Ey[opp]))]
    r = Ey[j] / Ey[i_peak]
    print(f"reflected peak: t={t[j]:.4e} s  Ey={Ey[j]:.6e} V/m")
    print(f"reflection coefficient r = {r:.6f}  |r| = {abs(r):.6f}  (expect r = -1 for PEC)")
else:
    print("no reflected pulse detected")

ts, Eys = load("diags/reducedfiles/fp_shadow.txt")
print(f"shadow region max|Ey| = {np.max(np.abs(Eys)):.6e} (should be ~0, PEC blocks field)")
