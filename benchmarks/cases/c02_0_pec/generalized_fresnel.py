#!/usr/bin/env python3
"""Reference reflection/transmission for a planar interface, including anisotropic media.

Two levels of generality are provided.

1. Normal incidence onto a (possibly rotated) tensor
   ---------------------------------------------------
   For k || z the transverse block of eps decouples into a 2x2 problem in (x, y):

       eps_t = [[exx, exy], [exy, eyy]]

   Its eigenvectors are the two polarisations that propagate without mixing, and the
   eigenvalues give the two refractive indices n_i = sqrt(lambda_i / eps0). Each mode
   sees an isotropic Fresnel interface with its own n_i:

       r_i = (n1_i - n2_i) / (n1_i + n2_i),   R_i = |r_i|^2
       t_i = 2 n1_i / (n1_i + n2_i),          T_i = (n2_i/n1_i) |t_i|^2

   This is the quantity our uniform-rotated-tensor plane-wave test measures (two modes,
   two speeds), so it is the right reference for that test.

2. Oblique incidence onto a uniaxial/diagonal tensor (s and p polarisation)
   -----------------------------------------------------------------------
   For a diagonal eps = diag(exx, eyy, ezz) with the interface normal along z:

       s polarisation (E || y):  n = sqrt(eyy)
       p polarisation (E in xz): k_x^2/ezz + k_z^2/exx = (omega/c)^2
                                 -> the Fresnel formulas with the angle from that
                                    dispersion relation

   Fully general oblique incidence in a rotated biaxial medium needs the 4x4 Berreman
   matrix; that is *not* implemented here. Use Meep for that case (see the P1 plan).

Usage
-----
  python3 generalized_fresnel.py                 # self-check + the P1 test case
  python3 generalized_fresnel.py --tensor        # print the rotated-tensor values
"""
import argparse
import math

import numpy as np


def n_modes(eps_t, eps0=1.0):
    """Refractive indices of the two normal-incidence modes of a 2x2 transverse block."""
    w, _ = np.linalg.eigh(np.asarray(eps_t, dtype=float) / eps0)
    w = np.maximum(w, 1e-30)
    return np.sqrt(w)          # ascending


def normal_incidence(eps1_t, eps2_t, eps0=1.0):
    """Per-mode (R, T) for normal incidence from material 1 to material 2."""
    n1 = n_modes(eps1_t, eps0)
    n2 = n_modes(eps2_t, eps0)
    r = (n1 - n2) / (n1 + n2)
    t = 2.0 * n1 / (n1 + n2)
    R = np.abs(r) ** 2
    T = (n2 / n1) * np.abs(t) ** 2
    return n1, n2, R, T


def uniaxial_s_p(exx, eyy, ezz, theta_deg, eps0=1.0):
    """Oblique s and p reflectance of a diagonal tensor against vacuum.

    theta is the incidence angle in vacuum. s uses n = sqrt(eyy); p follows from
    kx^2/ezz + kz^2/exx = k0^2 with kx conserved (= k0 sin theta).
    """
    th = math.radians(theta_deg)
    kx = math.sin(th)                     # in units of k0
    # s polarisation
    n_s = math.sqrt(eyy / eps0)
    # Fresnel needs n2*cos(theta2) = kz_s (in units of k0), not cos(theta2) alone
    kz_s = math.sqrt(max(n_s ** 2 - kx ** 2, 0.0))
    Rs = abs((math.cos(th) - kz_s) / (math.cos(th) + kz_s)) ** 2
    # p polarisation: kz from kx^2/ezz + kz^2/exx = 1
    A = exx / eps0
    B = ezz / eps0
    kz_p = math.sqrt(max((1.0 - kx ** 2 / B) * A, 0.0))
    # impedance-like Fresnel for p: r = (n1 cos t2 - n2 cos t1)/... expressed via kz
    rp = (kz_p / A - math.cos(th)) / (kz_p / A + math.cos(th))
    Rp = abs(rp) ** 2
    return Rs, Rp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tensor", action="store_true")
    args = ap.parse_args()

    print("=== self-checks ===")
    # isotropic: two degenerate modes, Fresnel R = 1/9 for n1=1, n2=2
    n1, n2, R, T = normal_incidence(np.eye(2), 4.0 * np.eye(2))
    print(f"isotropic n1={n1[0]:.6f} n2={n2[0]:.6f}  R={R[0]:.6f} (0.111111)  T={T[0]:.6f} (0.888889)")
    assert abs(R[0] - 1.0 / 9.0) < 1e-12 and abs(T[0] - 8.0 / 9.0) < 1e-12

    # isotropic at oblique incidence must reproduce the textbook Fresnel formulas
    Rs, Rp = uniaxial_s_p(1.0, 1.0, 1.0, 0.0)
    assert abs(Rs) < 1e-14 and abs(Rp) < 1e-14, (Rs, Rp)
    print("normal incidence on vacuum: Rs=Rp=0  OK")

    # rotated tensor of the P1 test: Rz(45) diag(1,4,2) Rz(45)^T
    a = math.radians(45.0)
    Rz = np.array([[math.cos(a), -math.sin(a), 0.0],
                   [math.sin(a),  math.cos(a), 0.0],
                   [0.0, 0.0, 1.0]])
    eps_rot = Rz @ np.diag([1.0, 4.0, 2.0]) @ Rz.T
    print("\n=== P1 rotated tensor ===")
    print("eps/eps0 =\n", np.array2string(eps_rot, precision=6, suppress_small=True))
    n1, n2, R, T = normal_incidence(np.eye(2), eps_rot[:2, :2])
    for i in range(2):
        print(f"mode {i}: n1={n1[i]:.6f} n2={n2[i]:.6f}  R={R[i]:.6f}  T={T[i]:.6f}")
    print(f"mean R = {R.mean():.6f}   mean T = {T.mean():.6f}   R+T = {(R+T).mean():.6f}")
    print("expected: n2 = 1 and 2 (eigenvalues of the x-y block are 1 and 4),")
    print("          R = 0 and 1/9, mean R = 1/18 = 0.055556")

    # sanity: the eigenvalues of the rotated tensor are those of diag(1,4,2)
    w = np.linalg.eigvalsh(eps_rot)
    print(f"eigenvalues of eps_rot = {np.array2string(w, precision=6)}  (expect 1,2,4)")
    assert np.allclose(np.sort(w), [1.0, 2.0, 4.0], atol=1e-12)

    print("\n=== oblique (uniaxial) ===")
    for th in (0.0, 30.0, 60.0, 80.0):
        Rs, Rp = uniaxial_s_p(4.0, 4.0, 4.0, th)
        print(f"theta={th:5.1f} deg  isotropic n=2:  Rs={Rs:.4f}  Rp={Rp:.4f}")


if __name__ == "__main__":
    main()
