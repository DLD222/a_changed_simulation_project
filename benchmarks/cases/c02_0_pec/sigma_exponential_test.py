#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""P2 reference: matrix-exponential update of a conductive anisotropic medium.

The semi-discrete field equation for a lossy dielectric is

    dE/dt = -A E + eps^-1 curlH ,      A = eps^-1 sigma          (1)

with curlH held constant over one step, the exact solution is

    E^{n+1} = exp(-A dt) E^n + [(1 - exp(-A dt)) / A] eps^-1 curlH     (2)

A is generally not symmetric, but with eps and sigma both SPD,

    B = eps^{-1/2} sigma eps^{-1/2}        (symmetric, positive semi-definite)
    A = eps^{-1/2} B eps^{1/2}

so both matrix functions reduce to functions of the symmetric B:

    exp(-A dt)          = eps^{-1/2} exp(-B dt) eps^{1/2}
    (1-exp(-A dt))/A    = eps^{-1/2} [(1-exp(-B dt))/B] eps^{1/2}

and f(B) = V diag(f(lambda_i)) V^T for B = V diag(lambda) V^T.

This script is the reference for the C++ port; it self-checks the limiting cases.

Usage:
  python3 sigma_exponential_test.py
"""
import numpy as np

EPS0 = 8.8541878128e-12


def sym(a):
    return 0.5 * (a + a.T)


def sqrtm_spd(a):
    w, v = np.linalg.eigh(sym(a))
    assert w.min() > 0, "matrix must be positive definite"
    return v @ np.diag(np.sqrt(w)) @ v.T


def invsqrtm_spd(a):
    w, v = np.linalg.eigh(sym(a))
    assert w.min() > 0, "matrix must be positive definite"
    return v @ np.diag(1.0 / np.sqrt(w)) @ v.T


def f_of_sym(b, kind, dt):
    """f(B) for a symmetric B; `kind` selects exp(-dt B) or (1-exp(-dt B))/B."""
    w, v = np.linalg.eigh(sym(b))
    if kind == "exp":
        g = np.exp(-w * dt)
    elif kind == "int":                     # (1 - e^{-lambda dt}) / lambda, -> dt as lambda->0
        g = np.where(np.abs(w) > 1e-300, (1.0 - np.exp(-w * dt)) / np.where(w == 0, 1.0, w), dt)
    else:
        raise ValueError(kind)
    return v @ np.diag(g) @ v.T


def exponential_update(E, curlH, eps, sigma, dt):
    """One step of eq. (2): E, curlH are 3-vectors; eps, sigma are 3x3 SPD."""
    eps_h = sqrtm_spd(eps)
    eps_ih = invsqrtm_spd(eps)
    B = sym(eps_ih @ sigma @ eps_ih)
    Edec = f_of_sym(B, "exp", dt)
    Eint = f_of_sym(B, "int", dt)
    S = eps_ih @ Eint @ eps_ih                      # multiplies curlH directly
    return Edec @ (eps_h @ E) * 0 + (eps_ih @ Edec @ eps_h) @ E + S @ curlH


def alpha_beta(sigma, eps, dt):
    """The scalar scheme currently in MacroscopicEvolveE (backward Euler)."""
    s = sigma * dt / eps
    return (1.0 - 0.5 * s) / (1.0 + 0.5 * s), dt / eps / (1.0 + 0.5 * s)


def main():
    dt = 1.0e-12
    I3 = np.eye(3)
    ok = True

    # --- [1] sigma = 0 must reduce to standard Yee: E + dt eps^-1 curlH
    eps = 4 * EPS0 * I3
    sigma = np.zeros((3, 3))
    E = np.array([1.0, 2.0, 3.0])
    curlH = np.array([0.5, -0.3, 0.8])
    En = exponential_update(E, curlH, eps, sigma, dt)
    ref = E + dt * np.linalg.inv(eps) @ curlH
    e1 = np.abs(En - ref).max()
    print(f"[1] sigma=0 -> Yee         max|diff| = {e1:.3e}")
    ok &= e1 < 1e-9 * max(1.0, np.abs(ref).max())

    # --- [2] scalar limit: compare with the backward-Euler alpha/beta scheme
    #     (they differ at O(dt^2), so the check is a small-dt convergence test)
    sig = 2.0            # tau = eps/sigma = 1.77e-11 s, so dt << tau below
    eps = 4 * EPS0 * I3
    sigma = sig * I3
    E = np.array([1.0, 1.0, 1.0])
    curlH = np.array([0.2, 0.2, 0.2])
    errs = []
    for dtk in (1e-12, 1e-13, 1e-14):   # dt/tau = 5.6e-2, 5.6e-3, 5.6e-4
        a, b = alpha_beta(sig, 4 * EPS0, dtk)
        Ebe = a * E + b * curlH
        Eex = exponential_update(E, curlH, eps, sigma, dtk)
        errs.append(np.abs(Ebe - Eex).max())
    ratio = errs[0] / errs[2]
    print(f"[2] scalar vs backward Euler  errs = {['%.2e' % e for e in errs]}  ratio(dt x100) = {ratio:.1f}")
    ok &= errs[0] < 1e-2 and errs[2] < errs[0] / 50

    # --- [3] homogeneous conductor relaxation: E(t) = E0 exp(-t/tau), tau = eps/sigma
    eps = 4 * EPS0 * I3
    sigma = 2.0e3 * I3
    tau = 4 * EPS0 / 2.0e3
    E = np.array([1.0, 0.0, 0.0])
    zero = np.zeros(3)
    nstep, dt = 400, tau / 400
    for _ in range(nstep):
        E = exponential_update(E, zero, eps, sigma, dt)
    exact = np.exp(-nstep * dt / tau)
    print(f"[3] relaxation tau={tau:.3e}s  E_num/E0={E[0]:.6f}  exact={exact:.6f}  rel={E[0]/exact-1:+.2e}")
    ok &= abs(E[0] / exact - 1) < 1e-10

    # --- [4] anisotropic sigma: each principal axis decays with its own tau
    R = np.array([[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]])
    sig_diag = np.diag([1.0e3, 4.0e3, 9.0e3])
    sigma = R @ sig_diag @ R.T
    eps = 4 * EPS0 * I3
    E = np.array([1.0, 0.0, 0.0])
    for _ in range(200):
        E = exponential_update(E, zero, eps, sigma, 4 * EPS0 / 9.0e3 / 200)
    # decompose along principal axes and compare each with exp(-sigma_i t/eps)
    t = 4 * EPS0 / 9.0e3
    amp0 = R.T @ np.array([1.0, 0.0, 0.0])
    amp = R.T @ E
    exact = amp0 * np.array([np.exp(-s * t / (4 * EPS0)) for s in np.diag(sig_diag)])
    print(f"[4] anisotropic decay       amp={np.round(amp,6)}  exact={np.round(exact,6)}")
    ok &= np.abs(amp - exact).max() < 1e-10

    # --- [5] energy monotone: E^T eps E must not increase without sources
    eps = 4 * EPS0 * I3
    sigma = np.array([[3e3, 1e3, 0], [1e3, 2e3, 0], [0, 0, 5e2]])
    E = np.array([1.0, -0.5, 0.3])
    prev = E @ eps @ E
    mono = True
    for _ in range(50):
        E = exponential_update(E, zero, eps, sigma, 1e-13)
        cur = E @ eps @ E
        if cur > prev * (1 + 1e-12):
            mono = False
        prev = cur
    print(f"[5] energy monotone         {'PASS' if mono else 'FAIL'}  final/initial={prev/(np.array([1.0,-0.5,0.3]) @ eps @ np.array([1.0,-0.5,0.3])):.4f}")
    ok &= mono

    # --- [6] stability at large sigma*dt: exp(-lambda dt) stays in (0,1]
    eps = 4 * EPS0 * I3
    sigma = 1e12 * I3
    for dtk in (1e-15, 1e-12, 1e-9):
        E = np.array([1.0, 1.0, 1.0])
        En = exponential_update(E, zero, eps, sigma, dtk)
        assert np.all(np.abs(En) <= np.abs(E) + 1e-12), "exponential must be bounded"
    print("[6] large sigma*dt bounded   PASS")

    print("\nRESULT:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
