#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C03-S5: lossy plane-wave attenuation from two point probes.

Theory (exact complex wavenumber, WarpX eps0):
    k = omega sqrt(mu eps) sqrt(1 - i sigma/(omega eps))
    |E(z2)/E(z1)| = exp(-|Im k| dz),   arg = -Re(k) dz  (mod 2 pi)

Usage:  python3 analyze.py diags/reducedfiles [--t0 20e-9] [--t1 95e-9]
"""
import argparse
import os
import numpy as np

EPS0 = 8.8541878128e-12
C = 299792458.0
MU0 = 1.0 / (EPS0 * C * C)


def load(path):
    d = np.loadtxt(path, skiprows=1)
    return d[:, 1], d[:, 6]          # time, Ey


def dft(t, y, w):
    return np.trapezoid(y * np.exp(-1j * w * t), t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diag", nargs="?", default="diags/reducedfiles")
    ap.add_argument("--eps-r", type=float, default=4.0)
    ap.add_argument("--sigma", type=float, default=0.01)
    ap.add_argument("--dz", type=float, default=0.5)
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--t1", type=float, default=None)
    a = ap.parse_args()

    t1, y1 = load(os.path.join(a.diag, "q1.txt"))
    t2, y2 = load(os.path.join(a.diag, "q2.txt"))
    tmax = t1.max() if a.t1 is None else a.t1
    m1 = (t1 >= a.t0) & (t1 <= tmax)
    m2 = (t2 >= a.t0) & (t2 <= tmax)

    # carrier from the probe spectrum
    y = y1[m1] - y1[m1].mean()
    f = np.fft.rfftfreq(len(y), t1[1] - t1[0])
    w = 2 * np.pi * f[np.argmax(np.abs(np.fft.rfft(y)))]
    print(f"FFT peak omega = {w:.6e}  (f = {w/(2*np.pi):.6e} Hz)")

    A1 = dft(t1[m1], y1[m1], w)
    A2 = dft(t2[m2], y2[m2], w)
    ratio = A2 / A1
    amp, ph = abs(ratio), np.angle(ratio) * 180 / np.pi

    eps = a.eps_r * EPS0
    k = w * np.sqrt(MU0 * eps) * np.sqrt(1 - 1j * a.sigma / (w * eps))
    kr, ki = k.real, abs(k.imag)
    amp_th = np.exp(-ki * a.dz)
    ph_th = np.angle(np.exp(-1j * kr * a.dz)) * 180 / np.pi
    print(f"measured |q2/q1| = {amp:.6f}   theory {amp_th:.6f}   "
          f"dev {100*(amp/amp_th-1):+.3f}%")
    print(f"measured phase   = {ph:+.3f} deg   theory {ph_th:+.3f} deg   "
          f"dev {ph-ph_th:+.3f} deg")
    print(f"k = {kr:.6f} {ki:+.6f}i   (1/m)")


if __name__ == "__main__":
    main()
