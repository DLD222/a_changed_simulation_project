#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Airy (Fabry-Perot) reference for a slab of thickness d, n2, between n1=n3.

  r12=(n1-n2)/(n1+n2), r23=(n2-n1)/(n2+n1)=-r12, t12=2n1/(n1+n2), t23=2n2/(n2+n1)
  delta = n2 * omega * d / c
  t = t12 t23 e^{i delta} / (1 + r12 r23 e^{2 i delta})
  R = |r|^2,  T = (n1/n1) |t|^2 = |t|^2

Self-checks: d=0 -> T=1; delta=m*pi -> T=1 (resonance); delta=m*pi+pi/2 -> T=1/(1+F).
Run: python3 airy_slab.py
"""
import numpy as np

C = 299792458.0

def airy(omega, n1, n2, d):
    r12 = (n1 - n2) / (n1 + n2)
    r23 = -r12
    t12 = 2.0 * n1 / (n1 + n2)
    t23 = 2.0 * n2 / (n2 + n1)
    delta = n2 * omega * d / C
    t = t12 * t23 * np.exp(1j * delta) / (1.0 + r12 * r23 * np.exp(2j * delta))
    r = (r12 + r23 * np.exp(2j * delta)) / (1.0 + r12 * r23 * np.exp(2j * delta))
    return np.abs(r) ** 2, np.abs(t) ** 2

if __name__ == "__main__":
    n1, n2 = 1.0, 2.0
    R12 = ((n1 - n2) / (n1 + n2)) ** 2
    F = 4.0 * R12 / (1.0 - R12) ** 2
    print("=" * 60)
    # d=0
    R, T = airy(2 * np.pi * C / 0.5, n1, n2, 0.0)
    print(f"[1] d=0: R={R:.6e} T={T:.6e} (expect 0,1)")
    # resonance delta = m*pi at omega = m*pi*c/(n2 d)
    d = 0.5
    for m in (1, 2, 3):
        w = m * np.pi * C / (n2 * d)
        R, T = airy(w, n1, n2, d)
        print(f"[2] resonance m={m}: R={R:.6e} T={T:.6e} (expect 1)")
    # anti-resonance delta = m*pi + pi/2
    w = (np.pi / 2) * C / (n2 * d)
    R, T = airy(w, n1, n2, d)
    print(f"[3] anti-resonance: T={T:.6e}  expect 1/(1+F)={1/(1+F):.6e}")
    # Fabry-Perot transmission spectrum peaks
    w = np.linspace(0.1, 2000.0, 200001) * 1e9
    R, T = airy(w, n1, n2, 1e-3)
    peaks = w[1:-1][(T[1:-1] > T[:-2]) & (T[1:-1] > T[2:])]
    print(f"[4] first three FP peaks at {peaks[:3]/1e9} GHz (analytic spacing "
          f"{np.pi*C/(n2*1e-3)/1e9:.4f} GHz)")
    print("=" * 60)
