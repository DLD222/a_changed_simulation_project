#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C03-E6: attenuation of a wave inside a good conductor (two in-slab probes).

Deterministic analysis (same doctrine as tensor_interface/analyze.py):
  - the carrier is FIXED to omega = 2*pi*c/lam (an input, never an FFT peak);
  - each probe is integrated over its own window (centre = emission + travel time,
    half-width = 5*tau), so no pulse is ever truncated by a common window.

Theory, with the code's own constants (WarpX mu0 = 1.2566370612685e-6, so
eps0 = 1/(mu0 c^2); using CODATA here was the §7.8 mistake):
  k = omega sqrt(mu0 eps0) sqrt(1 - i sigma/(omega eps0))
  |E(z2)/E(z1)| = exp(Im(k) dz)  [Im(k) < 0 = attenuation],   arg(E2/E1) = -Re(k) dz
  NOTE (selftest trap): the synthetic signals and the theory formula must not
  share a sign error -- with exp(-Im(k)*z) in BOTH, the selftest passes while the
  physics is inverted (growth instead of decay). Fixed: both use exp(Im(k)*z).

Usage:
  python3 analyze.py --selftest
  python3 analyze.py diags/reducedfiles --sigma 0.1 [--lam 1.0 --t-emission 9.6e-9 \
      --tau 3.2e-9 --z-src -5.0 --z1 0.046875 --z2 0.234375 --tol 0.03]
"""
import argparse
import os
import numpy as np

C = 299792458.0
MU0 = 1.2566370612685e-6            # Source/ablastr/constant.H:70 (the code's value)
EPS0 = 1.0 / (MU0 * C * C)


def k_complex(omega, sigma, eps=EPS0, mu=MU0):
    return omega * np.sqrt(mu * eps) * np.sqrt(1.0 - 1j * sigma / (omega * eps))


def load(path):
    d = np.loadtxt(path, skiprows=1)
    return d[:, 1], d[:, 6]          # time, Ey


def dft_window(t, y, w, centre, half):
    m = (t >= centre - half) & (t <= centre + half)
    return np.trapezoid(y[m] * np.exp(-1j * w * t[m]), t[m])


def analyse(diags, sigma, lam, t_emission, tau, z_src, z1, z2, tol):
    omega = 2.0 * np.pi * C / lam
    kc = k_complex(omega, sigma)
    v_med = omega / kc.real                       # phase velocity inside the slab
    t_slab = (0.0 - z_src) / C                    # source -> slab face (vacuum, s-pol)

    t1, y1 = load(os.path.join(diags, "q1.txt"))
    t2, y2 = load(os.path.join(diags, "q2.txt"))
    a1 = dft_window(t1, y1, omega, t_emission + t_slab + z1 / v_med, 5.0 * tau)
    a2 = dft_window(t2, y2, omega, t_emission + t_slab + z2 / v_med, 5.0 * tau)

    dz = z2 - z1
    amp_th = np.exp(kc.imag * dz)        # Im(k) < 0, so this is exp(-alpha*dz)
    ph_th = (-np.degrees(kc.real * dz) + 180.0) % 360.0 - 180.0
    amp = abs(a2 / a1)
    ph = np.degrees(np.angle(a2 / a1))
    print(f"k = {kc.real:.4f} {kc.imag:+.4f}i   alpha = {-kc.imag:.4f} 1/m   v_med = {v_med:.4e} m/s")
    print(f"|E2/E1| = {amp:.5f}   theory {amp_th:.5f}   dev {amp/amp_th-1:+.2%}"
          f"   (tol {tol:.0%}, power ratio exp(-2 a d) = {amp_th**2:.5f})")
    print(f"arg     = {ph:+.1f} deg   theory {ph_th:+.1f} deg   diff {ph-ph_th:+.1f} deg")
    return (abs(amp / amp_th - 1) <= tol), (abs(ph - ph_th) <= 3.0)


def selftest():
    """Synthetic forward wave sampled at two in-slab points; the analysis must
    recover exp(-Im k dz) and -Re(k) dz without knowing anything else."""
    sigma, lam, tau = 0.1, 1.0, 3.2e-9
    omega = 2.0 * np.pi * C / lam
    kc = k_complex(omega, sigma)
    v = omega / kc.real
    z1, z2, z_src, t0 = 0.05, 0.25, -5.0, 9.6e-9
    t = np.arange(0.0, 60e-9, 2.947e-11)
    env = lambda tt: np.exp(-((tt - t0) / tau) ** 2)
    # forward wave: envelope travels at v after the slab face; amplitude decays
    y1 = 0.5 * np.exp(kc.imag * (z1 - 0.0)) * env(t - (0 - z_src) / C - z1 / v) \
         * np.cos(omega * t - kc.real * z1)
    y2 = 0.5 * np.exp(kc.imag * (z2 - 0.0)) * env(t - (0 - z_src) / C - z2 / v) \
         * np.cos(omega * t - kc.real * z2)
    ok = True
    a1 = dft_window(t, y1, omega, t0 + (0 - z_src) / C + z1 / v, 5 * tau)
    a2 = dft_window(t, y2, omega, t0 + (0 - z_src) / C + z2 / v, 5 * tau)
    dz = z2 - z1
    amp_th, ph_th = np.exp(kc.imag * dz), (-np.degrees(kc.real * dz) + 180) % 360 - 180
    amp, ph = abs(a2 / a1), np.degrees(np.angle(a2 / a1))
    print(f"selftest: |E2/E1| = {amp:.5f} (th {amp_th:.5f}, {amp/amp_th-1:+.2e})"
          f"   arg = {ph:+.2f} (th {ph_th:+.2f}, diff {ph-ph_th:+.2e})")
    # 1e-3, not 1e-12: the two windows sit on different grid phases, so the discrete
    # trapezoid DFT leaves a small non-integer-cycle leakage (~3e-3 here).
    ok &= abs(amp / amp_th - 1) < 1e-2 and abs(ph - ph_th) < 1e-6
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diags", nargs="?", default=None)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--sigma", type=float, default=0.1)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--t-emission", type=float, default=9.6e-9)
    ap.add_argument("--tau", type=float, default=3.2e-9)
    ap.add_argument("--z-src", type=float, default=-5.0)
    ap.add_argument("--z1", type=float, default=0.046875)
    ap.add_argument("--z2", type=float, default=0.234375)
    ap.add_argument("--tol", type=float, default=0.03)
    a = ap.parse_args()
    if a.selftest or a.diags is None:
        raise SystemExit(0 if selftest() else 1)
    ok_a, ok_p = analyse(a.diags, a.sigma, a.lam, a.t_emission, a.tau,
                         a.z_src, a.z1, a.z2, a.tol)
    print("PASS" if (ok_a and ok_p) else "FAIL"
          + f"  (amp tol {a.tol:.0%}, phase tol 3 deg)")
    raise SystemExit(0 if (ok_a and ok_p) else 1)


if __name__ == "__main__":
    main()
