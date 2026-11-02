#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C03-E6 CW: two-way fit E(z) = A exp(ikz) + B exp(-ikz) from N point probes.

Why: a pulse + single-window DFT has ~1% systematic error here (the conductor is
strongly dispersive and the envelope deforms). The ln|E| slope also picks up the
lambda/2 standing wave produced by a ~0.5% numerical reflection, which biases the
fitted attenuation. A CW run + complex-amplitude DFT + two-way fit separates:

    Im(k)  -> attenuation   (1/delta, unbiased by the standing wave)
    |B/A|  -> numerical reflection at the sigma interface

Usage:
  python3 analyze_cw.py <diags_dir> --z-list 0.20 0.24 ... --sigma 0.1 --lam 1.0
                        [--t0 40e-9] [--t1 80e-9]
  python3 analyze_cw.py --selftest
"""
import argparse
import numpy as np

C = 299792458.0
EPS0 = 8.8541878188e-12
MU0 = 1.0 / (EPS0 * C * C)


def theory(sigma, lam):
    omega = 2.0 * np.pi * C / lam
    k = omega * np.sqrt(MU0 * EPS0) * np.sqrt(1.0 - 1j * sigma / (omega * EPS0))
    return k


def dft(t, y, w, t0, t1):
    m = (t >= t0) & (t <= t1)
    return np.sum(y[m] * np.exp(-1j * w * t[m])) * (t[1] - t[0])


def fixed_k_fit(z, e, k):
    """Given k, solve E = A exp(ikz) + B exp(-ikz) linearly -> unique A, B."""
    M = np.column_stack([np.exp(1j * k * z), np.exp(-1j * k * z)])
    sol, *_ = np.linalg.lstsq(M, e, rcond=None)
    A, B = sol
    model = M @ sol
    rms = np.sqrt(np.mean(np.abs(model - e) ** 2)) / max(np.mean(np.abs(e)), 1e-300)
    return A, B, rms


def two_way_fit(z, e, k0):
    """Least-squares fit of E(z) = A exp(ikz) + B exp(-ikz); A,B,k complex."""
    from scipy.optimize import least_squares

    def model(p):
        ar, ai, br, bi, kr, ki = p
        A, B, k = ar + 1j * ai, br + 1j * bi, kr + 1j * ki
        return A * np.exp(1j * k * z) + B * np.exp(-1j * k * z)

    def resid(p):
        d = model(p) - e
        return np.concatenate([d.real, d.imag])

    p0 = np.array([1.0, 0.0, 0.0, 0.0, k0.real, k0.imag])
    sol = least_squares(resid, p0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
    ar, ai, br, bi, kr, ki = sol.x
    A, B, k = ar + 1j * ai, br + 1j * bi, kr + 1j * ki
    rms = np.sqrt(np.mean(np.abs(model(sol.x) - e) ** 2)) / max(np.mean(np.abs(e)), 1e-300)
    return A, B, k, rms


def selftest():
    kt = theory(0.1, 1.0)
    z = np.array([0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56])
    A, B = 1.0 + 0.0j, 0.005 * np.exp(0.7j)      # |B/A| = 5e-3
    e = A * np.exp(1j * kt * z) + B * np.exp(-1j * kt * z)
    Af, Bf, kf, rms = two_way_fit(z, e, kt)
    print(f"selftest: k = {kf.real:.6f} {kf.imag:+.6f}i  (th {kt.real:.6f} {kt.imag:+.6f}i)")
    print(f"          |B/A| = {abs(Bf/Af):.6e} (th {abs(B/A):.6e})   rms {rms:.2e}")
    ok = (abs(kf - kt) / abs(kt) < 1e-10 and abs(abs(Bf / Af) / abs(B / A) - 1) < 1e-8)
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diags", nargs="?")
    ap.add_argument("--z-list", type=float, nargs="*", default=None)
    ap.add_argument("--sigma", type=float, default=0.1)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--t0", type=float, default=40e-9)
    ap.add_argument("--t1", type=float, default=80e-9)
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--fixed-k", action="store_true",
                    help="solve A,B linearly with k = k_theory (unique, no degeneracy)")
    a = ap.parse_args()

    if a.selftest:
        return selftest()

    kt = theory(a.sigma, a.lam)
    w = 2.0 * np.pi * C / a.lam
    z = np.asarray(a.z_list if a.z_list else
                   [0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56])
    e = []
    for i in range(len(z)):
        d = np.loadtxt(f"{a.diags}/p{i}.txt", skiprows=1)
        e.append(dft(d[:, 1], d[:, 6], w, a.t0, a.t1))
    e = np.asarray(e)

    if a.fixed_k:
        A, B, rms = fixed_k_fit(z, e, kt)
        k = kt
        print(f"fixed-k fit (k = k_theory):")
        print(f"  |A| = {abs(A):.6e}   |B| = {abs(B):.6e}   |B/A| = {abs(B/A):.6e}")
        print(f"  rms = {rms:.3e}")
        d = -np.log(abs(A) / abs(A)) if False else 0.0
        return 0 if rms < 0.05 else 1

    A, B, k, rms = two_way_fit(z, e, kt)
    print(f"theory k = {kt.real:.6f} {kt.imag:+.6f}i   1/delta = {-kt.imag:.6f}")
    print(f"fit    k = {k.real:.6f} {k.imag:+.6f}i   1/delta = {-k.imag:.6f}")
    print(f"  Im(k) rel dev = {(k.imag - kt.imag)/abs(kt.imag):+.3%}   "
          f"Re(k) rel dev = {(k.real - kt.real)/abs(kt.real):+.3%}")
    print(f"  |B/A| = {abs(B/A):.6e}   fit rms {rms:.2e}")
    ok = abs((k.imag - kt.imag) / abs(kt.imag)) <= a.tol
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
