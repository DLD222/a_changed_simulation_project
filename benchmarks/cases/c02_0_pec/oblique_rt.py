#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Oblique-incidence R/T from two probe planes per side.

Given complex amplitudes E1, E2 at z1, z2 (incident side, z1<z2<interface) and
E3, E4 at z3, z4 (transmitted side), plus the (measured) normal wavenumbers
k1 (incident side) and k2 (transmitted side), decompose

  incident side:  E(z) = Ei e^{i k1 z} + Er e^{-i k1 z}
  transmitted:    E(z) = Et e^{i k2 z}

and return R=|Er/Ei|^2, T=(k2/k1)|Et/Ei|^2 for s polarization (for p use the
H-field amplitudes or the full Poynting flux).

Self-test with synthetic Ei, Er, Et.
Run: python3 oblique_rt.py
"""
import numpy as np


def decompose(E1, E2, z1, z2, k1):
    """Return (Ei, Er) from two complex amplitudes on the incident side."""
    D = np.exp(1j * k1 * (z1 - z2)) - np.exp(1j * k1 * (z2 - z1))
    Ei = (E1 * np.exp(-1j * k1 * z2) - E2 * np.exp(-1j * k1 * z1)) / D
    Er = (E2 * np.exp(1j * k1 * z1) - E1 * np.exp(1j * k1 * z2)) / D
    return Ei, Er


def transmitted(E3, E4, z3, z4, k2):
    """Return (Et, consistency) from two planes on the transmitted side."""
    Et3 = E3 * np.exp(-1j * k2 * z3)
    Et4 = E4 * np.exp(-1j * k2 * z4)
    return 0.5 * (Et3 + Et4), abs(Et3 - Et4) / max(abs(Et3), 1e-300)


def measure_k(E1, E2, z1, z2):
    """k from a single forward wave at two planes (unwrapped phase)."""
    return np.angle(E2 * np.conj(E1)) / (z2 - z1)


if __name__ == "__main__":
    print("=" * 60)
    k1, k2 = 8.0, 12.0
    z1, z2, z3, z4 = -0.6, -0.3, 0.3, 0.6
    r = (k1 - k2) / (k1 + k2); t = 2.0 * k1 / (k1 + k2)
    Ei, Er, Et = 1.0 + 0.0j, r + 0.0j, t + 0.0j
    E1 = Ei * np.exp(1j * k1 * z1) + Er * np.exp(-1j * k1 * z1)
    E2 = Ei * np.exp(1j * k1 * z2) + Er * np.exp(-1j * k1 * z2)
    E3 = Et * np.exp(1j * k2 * z3)
    E4 = Et * np.exp(1j * k2 * z4)
    # k measurement from a pure forward run (Er=0)
    kk = measure_k(Ei * np.exp(1j * k1 * z1), Ei * np.exp(1j * k1 * z2), z1, z2)
    print(f"[1] measure_k: {kk:.6f} (expect {k1})")
    Eii, Err = decompose(E1, E2, z1, z2, k1)
    print(f"[2] decompose: Ei={Eii:.9f} Er={Err:.9f} (expect {Ei}, {Er})")
    Ett, cons = transmitted(E3, E4, z3, z4, k2)
    print(f"[3] transmitted: Et={Ett:.9f} consistency={cons:.3e} (expect {Et})")
    R = abs(Err / Eii) ** 2
    T = (k2 / k1) * abs(Ett / Eii) ** 2
    print(f"[4] R={R:.9f} T={T:.9f} R+T={R+T:.9f}")
    print("=" * 60)
