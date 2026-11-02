#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Evanescent 1D discrete dispersion by FORWARD mapping (no implicit root finding).

State (u = Ey^n at a node, v = Bx^{n-1/2} at a face). One step is a 2x2 matrix M(k):

    M11 = exp(-a dt) + C (1 - e^{-ik dz}) (dt/dz) (e^{ik dz} - 1)
    M12 = C (1 - e^{-ik dz})
    M21 = (dt/dz) (e^{ik dz} - 1)
    M22 = 1
    C   = (1 - exp(-a dt)) / (a mu eps dz),   a = sigma / eps

For a given complex k: take the eigenvalue of M closest to exp(-i omega0 dt), then
omega(k) = i ln(lambda) / dt.  The forward map k -> omega is single-valued near the
physical branch, so solving omega(k) = omega0 with bounded least_squares cannot fall
onto the 2x/10x spurious roots that fsolve hit earlier.

Usage:
  python3 dispersion_scan.py --compare
  python3 dispersion_scan.py --dz 0.03125 --dt 2.948317921e-11 --sigma 0.1
"""
import argparse
import numpy as np
from scipy.optimize import least_squares

C_LIGHT = 299792458.0
EPS0 = 8.8541878188e-12
MU0 = 1.0 / (EPS0 * C_LIGHT * C_LIGHT)


def omega_of_k(k, omega0, sigma, eps, mu, dz, dt):
    a = sigma / eps
    C = (1.0 - np.exp(-a * dt)) / (a * mu * eps * dz)
    eik = np.exp(1j * k * dz)
    M = np.array([
        [np.exp(-a * dt) + C * (1 - eik) * (dt / dz) * (eik - 1), C * (1 - eik)],
        [(dt / dz) * (eik - 1), 1.0 + 0j],
    ])
    lam = np.linalg.eigvals(M)
    target = np.exp(-1j * omega0 * dt)
    lam = lam[np.argmin(np.abs(lam - target))]
    return 1j * np.log(lam) / dt


def k_discrete(omega0, sigma, eps, dz, dt, mu=MU0):
    kc = omega0 * np.sqrt(mu * eps) * np.sqrt(1.0 - 1j * sigma / (omega0 * eps))

    def resid(x):
        k = x[0] + 1j * x[1]
        w = omega_of_k(k, omega0, sigma, eps, mu, dz, dt)
        # scale: omega is ~1e9, so normalise to make the two rows comparable
        return np.array([(w.real - omega0) / omega0, w.imag / omega0])

    lo = [0.8 * kc.real, 1.2 * kc.imag]          # kc.imag < 0
    hi = [1.2 * kc.real, 0.8 * kc.imag]
    sol = least_squares(resid, [kc.real, kc.imag], bounds=(lo, hi),
                        xtol=1e-15, ftol=1e-15, gtol=1e-15)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.x[0] + 1j * sol.x[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dz', type=float, default=0.03125)
    ap.add_argument('--dt', type=float, default=2.948317921e-11)
    ap.add_argument('--sigma', type=float, default=0.1)
    ap.add_argument('--eps', type=float, default=EPS0)
    ap.add_argument('--lam', type=float, default=1.0)
    ap.add_argument('--compare', action='store_true')
    a = ap.parse_args()

    omega0 = 2.0 * np.pi * C_LIGHT / a.lam
    kc = omega0 * np.sqrt(MU0 * a.eps) * np.sqrt(1.0 - 1j * a.sigma / (omega0 * a.eps))
    print(f'continuum: k = {kc.real:.8f} {kc.imag:+.8f}i   '
          f'1/delta = {-kc.imag:.8f}')

    grid = ((0.03125, 2.948317921e-11), (0.015625, 1.4741589605e-11))
    if not a.compare:
        grid = ((a.dz, a.dt),)

    for dz, dt in grid:
        kd = k_discrete(omega0, a.sigma, a.eps, dz, dt)
        print(f'  dz={dz:.6f}: k_num = {kd.real:.8f} {kd.imag:+.8f}i   '
              f'1/delta_num = {-kd.imag:.8f}   '
              f'dev = {(-kd.imag) / (-kc.imag) - 1:+.4%}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
