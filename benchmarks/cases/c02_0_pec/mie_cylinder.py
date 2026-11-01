#!/usr/bin/env python3
"""2D Mie series for a PEC circular cylinder, plane wave along +z.

Convention: time dependence e^{-j w t}, so outgoing waves are Hankel functions of
the first kind H_n^(1). With theta measured from +z (the incidence direction),

  E parallel to the axis (WarpX: Ey, "TM^z" in the cylinder literature)
      E_y^inc = E0 e^{jkz} = E0 sum_n j^n J_n(kr) e^{jn theta}
      E_y^sc = -E0 sum_n j^n a_n H_n^(1)(kr) e^{jn theta},
      a_n = J_n(ka) / H_n^(1)(ka)

  E perpendicular to the axis (WarpX: By, "TE^z")
      B_y^inc = B0 sum_n j^n J_n(kr) e^{jn theta}
      B_y^sc = -B0 sum_n j^n b_n H_n^(1)(kr) e^{jn theta},
      b_n = J_n'(ka) / H_n^(1)'(ka)

The 2D scattering width (RCS per unit length, metres) is

      sigma(phi) = (2 pi / k) | sum_n j^n c_n e^{jn phi} |^2 ,  c_n = a_n or b_n

The series is truncated with the standard rule n_max = ceil(ka + 4 (ka)^(1/3) + 8),
which is generous for PEC cylinders.

Run as a script for a self-check:
    python3 mie_cylinder.py
"""
import numpy as np
from scipy.special import hankel1, jv


def nmax(ka):
    return int(np.ceil(ka + 4.0 * ka ** (1.0 / 3.0) + 8.0))


def coeffs_tmz(ka, n_max=None):
    """a_n for E parallel to the axis (Ey)."""
    n = nmax(ka) if n_max is None else n_max
    ns = np.arange(-n, n + 1)
    return ns, jv(ns, ka) / hankel1(ns, ka)


def coeffs_tez(ka, n_max=None):
    """b_n for E perpendicular to the axis (By), derivative ratio."""
    n = nmax(ka) if n_max is None else n_max
    ns = np.arange(-n, n + 1)
    jp = 0.5 * (jv(ns - 1, ka) - jv(ns + 1, ka))
    hp = 0.5 * (hankel1(ns - 1, ka) - hankel1(ns + 1, ka))
    return ns, jp / hp


def field(k, a, r, theta, pol, n_max=None):
    """Complex amplitude of the total field (incident + scattered), unit input.

    pol = 'ey' (E parallel to the axis) or 'by' (E perpendicular, B parallel).
    Returns an array shaped like r/theta (broadcast).
    """
    ka = k * a
    ns, c = coeffs_tmz(ka, n_max) if pol == "ey" else coeffs_tez(ka, n_max)
    r = np.asarray(r, dtype=float)
    theta = np.asarray(theta, dtype=float)
    out = np.zeros(np.broadcast(r, theta).shape, dtype=complex)
    for i, n in enumerate(ns):
        phase = np.exp(1j * n * theta)
        out += (1j ** n) * (jv(n, k * r) - c[i] * hankel1(n, k * r)) * phase
    return out


def scattered(k, a, r, theta, pol, n_max=None):
    """Complex amplitude of the scattered field only (unit incident amplitude)."""
    ka = k * a
    ns, c = coeffs_tmz(ka, n_max) if pol == "ey" else coeffs_tez(ka, n_max)
    r = np.asarray(r, dtype=float)
    theta = np.asarray(theta, dtype=float)
    out = np.zeros(np.broadcast(r, theta).shape, dtype=complex)
    for i, n in enumerate(ns):
        out -= (1j ** n) * c[i] * hankel1(n, k * r) * np.exp(1j * n * theta)
    return out


def rcs(k, a, phi, pol, n_max=None):
    """2D scattering width sigma(phi) in metres, normalised to unit incidence."""
    ka = k * a
    ns, c = coeffs_tmz(ka, n_max) if pol == "ey" else coeffs_tez(ka, n_max)
    phi = np.asarray(phi, dtype=float)
    s = np.zeros(phi.shape, dtype=complex)
    for i, n in enumerate(ns):
        s += (1j ** n) * c[i] * np.exp(1j * n * phi)
    # sigma = lim 2*pi*r*|E_sc|^2 with |E_i| = 1; the far-field form is
    # (4/k)|sum c_n e^{jn phi}|^2 (verified numerically in __main__)
    return (4.0 / k) * np.abs(s) ** 2


if __name__ == "__main__":
    # self-checks
    k = 2.0 * np.pi            # lambda = 1 m
    a = 1.0
    ka = k * a
    ns, An = coeffs_tmz(ka)
    _, Bn = coeffs_tez(ka)
    print(f"ka = {ka:.6f}, n_max = {nmax(ka)}, series terms = {len(ns)}")
    print(f"max|a_n| = {np.max(np.abs(An)):.6e}, max|b_n| = {np.max(np.abs(Bn)):.6e}")

    # PEC boundary: total field must vanish on r = a
    th = np.linspace(0, 2 * np.pi, 361)
    for pol in ("ey", "by"):
        tot = field(k, a, np.ones_like(th) * (a * (1 + 1e-9)), th, pol)
        print(f"pol={pol}: max|total(r=a)| = {np.max(np.abs(tot)):.3e} "
              f"(should be ~0 at the PEC surface)")

    # energy: <sigma> over phi must equal (4/k) sum |a_n|^2
    phi = np.linspace(0, 2 * np.pi, 2001)
    sig = rcs(k, a, phi, "ey")
    integ = np.trapezoid(sig, phi) / (2 * np.pi)
    total = (4.0 / k) * np.sum(np.abs(An) ** 2)
    print(f"<sigma> over phi = {integ:.6e} m, (4/k) sum|a_n|^2 = {total:.6e} m, "
          f"rel.diff = {abs(integ - total) / total:.3e}")

    # far-field consistency: 2 pi r |E_sc|^2 must approach sigma(phi)
    r_far = 1.0e3
    for pol, c in (("ey", An), ("by", Bn)):
        sc = scattered(k, a, np.full_like(phi, r_far), phi, pol)
        lhs = 2 * np.pi * r_far * np.abs(sc) ** 2
        rhs = rcs(k, a, phi, pol)
        print(f"pol={pol}: max relative difference far field vs sigma = "
              f"{np.max(np.abs(lhs - rhs)) / np.max(rhs):.3e}")

    # TE^z boundary condition is dH_z/dr = 0 on the PEC surface (not H_z = 0)
    dr = 1e-6 * a
    h1 = field(k, a, np.full_like(th, a + dr), th, "by")
    h0 = field(k, a, np.full_like(th, a - dr), th, "by")
    print(f"pol=by: max|dH/dr| at r=a = {np.max(np.abs(h1 - h0)) / (2 * dr):.3e}")
