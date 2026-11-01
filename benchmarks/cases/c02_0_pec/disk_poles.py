#!/usr/bin/env python3
"""T7 (analytic reference): complex poles of the TM Mie denominator for a dielectric disk.

For a cylinder of radius a, index m, the TM (E || axis) coefficients have the
denominator

    D_l(x) = m J_l'(m x) H_l(x) - H_l'(x) J_l(m x),      x = k a

A resonance is a complex root D_l(x) = 0 with Im(x) < 0 (decaying mode, e^{-i w t}
convention). Its frequency and quality factor are

    f = Re(x) / (2 pi a),      Q = Re(x) / (2 |Im(x)|)

Usage:
  python3 disk_poles.py --m 3.4 --a 1.0 --l 3 4 5 --guess 5.9 7.4 8.8
"""
import argparse

import numpy as np
from scipy.optimize import newton
from scipy.special import hankel1, jv


def denom(x, m, l):
    J = jv(l, x)
    Jp = 0.5 * (jv(l - 1, x) - jv(l + 1, x))
    Jm = jv(l, m * x)
    Jpm = 0.5 * (jv(l - 1, m * x) - jv(l + 1, m * x))
    H = hankel1(l, x)
    Hp = 0.5 * (hankel1(l - 1, x) - hankel1(l + 1, x))
    return m * Jpm * H - Hp * Jm


def pole(m, l, guess, a=1.0):
    """Newton on the complex plane, with a numerical derivative (scipy handles it)."""
    root = newton(lambda x: denom(x, m, l), guess, tol=1e-12, maxiter=200)
    f = root.real / (2 * np.pi * a)
    Q = root.real / (2 * abs(root.imag))
    return root, f, Q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=float, default=3.4)
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--l", type=int, nargs="+", default=[3, 4, 5])
    ap.add_argument("--guess", type=float, nargs="+", default=[1.47, 1.84, 2.20])
    args = ap.parse_args()

    print(f"m={args.m}, a={args.a}")
    print(" l        x (complex)          f          Q")
    for l, g in zip(args.l, args.guess):
        try:
            root, f, Q = pole(args.m, l, complex(g, -0.02), args.a)
            print(f" {l}   {root.real:9.5f} {root.imag:+9.5f}i   {f:.5f}   {Q:8.1f}")
        except Exception as exc:     # noqa: BLE001
            print(f" {l}   failed: {exc}")


if __name__ == "__main__":
    main()
