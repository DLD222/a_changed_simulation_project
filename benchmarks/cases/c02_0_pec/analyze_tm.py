#!/usr/bin/env python3
"""TM normal-incidence PEC reflection coefficient from FieldProbe point data.
columns: 0 step, 1 time, 2 x, 3 y, 4 z, 5 Ex, 6 Ey, 7 Ez, 8 Bx, 9 By, 10 Bz, 11 S
Incident pulse arrives first (positive Ex), reflected pulse later (negative Ex, PEC r=-1)."""
import numpy as np, sys

path = sys.argv[1] if len(sys.argv) > 1 else "diags/reducedfiles/fp_inc.txt"
d = np.loadtxt(path, skiprows=1)
t = d[:, 1]
Ex = d[:, 5]

ip = int(np.argmax(Ex))                    # incident peak (positive)
t_inc = t[ip]
# search the reflected (negative) peak in a time window after the incident
win = (t > t_inc + 4e-9) & (t < t_inc + 60e-9)
if win.any():
    idx = np.where(win)[0]
    ineg = idx[np.argmin(Ex[win])]
    r = Ex[ineg] / Ex[ip]
    print(f"incident : t={t_inc*1e9:.3f} ns  Ex={Ex[ip]:+.6e} V/m")
    print(f"reflected: t={t[ineg]*1e9:.3f} ns  Ex={Ex[ineg]:+.6e} V/m")
    print(f"r = {r:+.6f}   |r| = {abs(r):.6f}   (PEC expects -1)")
else:
    print(f"incident : t={t_inc*1e9:.3f} ns  Ex={Ex[ip]:+.6e} V/m")
    print("no reflected pulse in window")
