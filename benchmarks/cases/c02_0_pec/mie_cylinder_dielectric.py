#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""2D Mie series for a *dielectric* circular cylinder, plane wave along +z.

Same convention as mie_cylinder.py (time dependence e^{-j w t}, H_n^(1) outgoing,
theta measured from +z):

  E || axis (WarpX Ey, "TM^z"):
      a_n = [ m J_n'(m x) J_n(x) - J_n'(x) J_n(m x) ]
            / [ m J_n'(m x) H_n(x) - H_n'(x) J_n(m x) ]
  H || axis (WarpX By, "TE^z"):
      b_n = [ J_n'(m x) J_n(x) - m J_n'(x) J_n(m x) ]
            / [ J_n'(m x) H_n(x) - m H_n'(x) J_n(m x) ]

with x = k a, m = sqrt(eps_r) (mu = mu0).  m -> 1 gives c_n = 0 (vacuum);
for a lossless dielectric the optical theorem holds: |c_n|^2 = Re(c_n),
so sigma_sca = sigma_ext = (4/k) sum |c_n|^2.

Self-checks in __main__:
  1. m=1 -> coefficients vanish
  2. optical theorem  |c_n|^2 = Re(c_n)   (lossless)
  3. large-x extinction paradox  sigma_ext/(4a) -> 1
Run:  python3 mie_cylinder_dielectric.py
"""
import numpy as np
from scipy.special import hankel1, jv


def nmax(x):
    return int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 8.0))


def _jp(ns, x):
    return 0.5 * (jv(ns - 1, x) - jv(ns + 1, x))


def _hp(ns, x):
    return 0.5 * (hankel1(ns - 1, x) - hankel1(ns + 1, x))


def coeffs_tmz(x, m, n_max=None):
    n = nmax(x) if n_max is None else n_max
    ns = np.arange(-n, n + 1)
    J, Jp = jv(ns, x), _jp(ns, x)
    Jm, Jpm = jv(ns, m * x), _jp(ns, m * x)
    H, Hp = hankel1(ns, x), _hp(ns, x)
    return ns, (m * Jpm * J - Jp * Jm) / (m * Jpm * H - Hp * Jm)


def coeffs_tez(x, m, n_max=None):
    n = nmax(x) if n_max is None else n_max
    ns = np.arange(-n, n + 1)
    J, Jp = jv(ns, x), _jp(ns, x)
    Jm, Jpm = jv(ns, m * x), _jp(ns, m * x)
    H, Hp = hankel1(ns, x), _hp(ns, x)
    return ns, (Jpm * J - m * Jp * Jm) / (Jpm * H - m * Hp * Jm)


def coeffs(x, m, pol, n_max=None):
    return coeffs_tmz(x, m, n_max) if pol == "ey" else coeffs_tez(x, m, n_max)


def scattered(x, m, kr, theta, pol, n_max=None):
    """Scattered complex amplitude, unit incident amplitude."""
    ns, c = coeffs(x, m, pol, n_max)
    kr = np.asarray(kr, float); theta = np.asarray(theta, float)
    out = np.zeros(np.broadcast(kr, theta).shape, complex)
    for i, n in enumerate(ns):
        out -= (1j ** n) * c[i] * hankel1(n, kr) * np.exp(1j * n * theta)
    return out


def rcs(x, m, k, phi, pol, n_max=None):
    """2D scattering width sigma(phi) [m], unit incidence."""
    ns, c = coeffs(x, m, pol, n_max)
    phi = np.asarray(phi, float)
    s = np.zeros(phi.shape, complex)
    for i, n in enumerate(ns):
        s += (1j ** n) * c[i] * np.exp(1j * n * phi)
    return (4.0 / k) * np.abs(s) ** 2


if __name__ == "__main__":
    print("=" * 72)
    # 1. m=1 -> vacuum, coefficients must vanish
    for pol in ("ey", "by"):
        x = 2.0
        ns, c = coeffs(x, 1.0, pol)
        print(f"[1] m=1, pol={pol}: max|c_n| = {np.max(np.abs(c)):.3e}")

    # 2. optical theorem for a lossless dielectric: |c_n|^2 = Re(c_n)
    for m in (1.5, 2.0, 4.0):
        for x in (0.5, 2.0, 5.0):
            for pol in ("ey", "by"):
                ns, c = coeffs(x, m, pol)
                err = np.max(np.abs(np.abs(c) ** 2 - np.real(c)))
                print(f"[2] m={m}, x={x}, pol={pol}: "
                      f"max||c|^2-Re(c)| = {err:.3e}")

    # 3. large-x extinction paradox: sigma_ext/(4a) -> 1
    for pol in ("ey", "by"):
        x = 100.0
        ns, c = coeffs(x, 2.0, pol)
        sigma = (4.0 / x) * np.sum(np.abs(c) ** 2) * x  # (4/k) sum|c|^2 with k=x/a
        # (4/k) sum |c|^2 = (4a/x) sum |c|^2
        sig = (4.0 / x) * np.sum(np.abs(c) ** 2)
        print(f"[3] x=100, pol={pol}: sigma_sca/(4a) = {sig/4.0:.4f}")
    print("=" * 72)
