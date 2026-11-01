#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Independent reference checks for the medium-solver math layer.

This exercises the formulas that `SubpixelTensor.H` implements, using numpy/scipy/
mpmath as the oracle.  It does not need a WarpX build.

  M1 isotropic chi^-1        -> harmonic (transverse) / arithmetic (normal) limits
  M2 anisotropic tau         -> same limits for a diagonal tensor
  M3 general 3D frame        -> orthonormal with e_z = n
  M4 expm M1/M2              -> scipy.linalg.expm (tensor eps, tensor sigma)
  M5 (1-e^-x)/x series       -> mpmath at 50 digits
"""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]          # benchmarks/cases
sys.path.insert(0, str(ROOT / "c02_0_pec"))
import p1_tensor_reference as ref                   # noqa: E402

EPS0 = 8.8541878128e-12

fails = 0


def check(name, got, want, tol):
    global fails
    got = np.asarray(got, dtype=float)
    want = np.asarray(want, dtype=float)
    err = float(np.max(np.abs(got - want)))
    ok = err <= tol
    fails += 0 if ok else 1
    print(f"{'PASS' if ok else 'FAIL'}  {name:42s} max|diff|={err:.3e}  (tol {tol:.0e})")


def chi_iso(eps1, eps2, f, n):
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    ea = f * eps1 + (1 - f) * eps2
    ei = f / eps1 + (1 - f) / eps2
    return np.outer(n, n) * (ei - 1 / ea) + np.eye(3) / ea


def expm_mats(eps, sigma, dt):
    """Same route as SubpixelTensor.H: eigh(eps) -> B = eps^-1/2 sigma eps^-1/2 ->
    eigh(B) -> M1 = eps^-1/2 exp(-B dt) eps^1/2, M2 = eps^-1/2 ((1-exp(-B dt))/B) eps^-1/2."""
    lam, V = np.linalg.eigh(eps)
    def fun(g):
        return V @ np.diag(g(lam)) @ V.T
    eps_h = fun(np.sqrt)
    eps_ih = fun(lambda l: 1 / np.sqrt(l))
    B = eps_ih @ sigma @ eps_ih
    lb, Vb = np.linalg.eigh(B)
    def fb(g):
        return Vb @ np.diag(g(lb)) @ Vb.T
    edec = fb(lambda l: np.exp(-l * dt))
    # (1 - e^{-x})/x with the series for small x
    def eint_f(l):
        x = l * dt
        return np.where(np.abs(x) < 1e-8, 1 - x / 2 + x * x / 6,
                        (1 - np.exp(-x)) / np.where(x == 0, 1, x))
    eint = dt * fb(eint_f)   # (1-e^{-lam dt})/lam = dt * (1-e^{-x})/x
    return eps_ih @ edec @ eps_h, eps_ih @ eint @ eps_ih


def main():
    rng = np.random.default_rng(12345)

    # ---- M1: isotropic chi^-1 limits, normal along z --------------------------------
    e1, e2, f = 4 * EPS0, EPS0, 0.37
    chi = chi_iso(e1, e2, f, [0, 0, 1])
    check("M1 chi_zz = harmonic mean (normal)", chi[2, 2], f / e1 + (1 - f) / e2, 1e-14)
    check("M1 chi_yy = 1/arithmetic (transverse)", chi[1, 1], 1 / (f * e1 + (1 - f) * e2), 1e-14)
    check("M1 chi_xy = 0", chi[0, 1], 0.0, 1e-15)

    # ---- M2: anisotropic tau -------------------------------------------------------
    # (a) spherical tensors must reduce to the isotropic closed form exactly
    T1iso = np.diag([4.0, 4.0, 4.0]) * EPS0
    T2iso = np.diag([1.0, 1.0, 1.0]) * EPS0
    chi_a = ref.chi_inv_anisotropic(T1iso, T2iso, f, [0, 0, 1])
    check("M2 spherical tensor = isotropic", chi_a, chi_iso(4 * EPS0, EPS0, f, [0, 0, 1]), 1e-13)
    # (b) diagonal tensor, n = z: eigenvalues = {harmonic(e_zz), 1/arith(e_xx), 1/arith(e_yy)}
    T1 = np.diag([4.0, 5.0, 6.0]) * EPS0
    T2 = np.diag([1.0, 1.0, 1.0]) * EPS0
    chi_d = ref.chi_inv_anisotropic(T1, T2, f, [0, 0, 1])
    ev = np.sort(np.linalg.eigvalsh(chi_d))
    want = np.sort([f / T1[2, 2] + (1 - f) / T2[2, 2],
                    1 / (f * T1[0, 0] + (1 - f) * T2[0, 0]),
                    1 / (f * T1[1, 1] + (1 - f) * T2[1, 1])])
    check("M2 diag tensor eigenvalues", ev, want, 1e-12)

    # ---- M3: general 3D frame is orthonormal with e_z = n --------------------------
    for k in range(5):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        R = ref.frame(n)
        check(f"M3 frame orthonormal {k}", R.T @ R, np.eye(3), 1e-14)
        check(f"M3 frame first column = n {k}", R[:, 0], n, 1e-14)

    # ---- M4: expm M1/M2 vs scipy for a random SPD eps and symmetric sigma ----------
    for k in range(5):
        eps = ref.rand_spd()
        sigma = ref.rand_spd() * 0.05                   # SPD, full rank
        dt = 1.0e-3
        M1, M2 = expm_mats(eps, sigma, dt)
        A = np.linalg.solve(eps, sigma)
        M1_ref = expm(-A * dt)
        M2_ref = (np.eye(3) - expm(-A * dt)) @ np.linalg.inv(A) @ np.linalg.inv(eps)
        check(f"M4 M1 = exp(-A dt) {k}", M1, M1_ref, 1e-13)
        check(f"M4 M2 = (1-exp(-A dt)) A^-1 eps^-1 {k}", M2, M2_ref, 1e-13)

    # ---- M5: (1-e^-x)/x series at small x -----------------------------------------
    mp.mp.dps = 50
    for x in (1e-8, 1e-10, 1e-12, 1e-15):
        exact = mp.mpf(1) - mp.e ** (-mp.mpf(x))
        exact = exact / mp.mpf(x)
        naive = (1.0 - np.exp(-x)) / x
        series = 1 - x / 2 + x * x / 6
        e_n = abs(naive - float(exact)) / float(exact)
        e_s = abs(series - float(exact)) / float(exact)
        ok = e_s <= 1e-15
        global fails
        fails += 0 if ok else 1
        print(f"{'PASS' if ok else 'FAIL'}  M5 series x={x:.0e}: naive {e_n:.1e}, "
              f"series {e_s:.1e}")

    print()
    print("ALL PASS" if fails == 0 else f"{fails} FAILURE(S)")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
