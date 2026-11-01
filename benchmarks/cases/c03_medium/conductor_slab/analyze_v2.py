#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""C03-E6 V2: skin-depth slope fit from N point probes inside the conductor.

Theory (good conductor, WarpX eps0):
    k = omega sqrt(mu eps) sqrt(1 - i sigma/(omega eps))     (Im k < 0)
    |E(z)| = |E0| exp(z * Im k) = |E0| exp(-z/delta),  delta = -1/Im(k)

Method (deterministic, same spirit as tensor_interface/analyze.py):
  * fixed carrier omega = 2 pi c / lam   (never FFT-peak; it drifts on truncation)
  * one DFT per probe over the SAME transmission window
  * fit ln|E(z)| = ln|E0| - z/delta over points with z >= 1.5 delta
    (the first ~1.5 delta is contaminated by the interface transient)

Usage:
  python3 analyze_v2.py <diags_dir> [--sigma 0.1] [--lam 1.0]
                        [--t0 30e-9] [--t1 55e-9] [--tol 0.03]
  python3 analyze_v2.py --selftest
"""
import argparse
import numpy as np

C = 299792458.0
EPS0 = 8.8541878188e-12          # WarpX ablastr::constant::epsilon_0_v (CODATA 2018)
MU0 = 1.0 / (EPS0 * C * C)

# probe z positions are CELL CENTRES of the E6 grid: z = zlo + (j + 0.5) * dz
# defaults are the 640^2 case (dz = 0.03125, j = 323,325,...,335)
Z_JS = (323, 325, 327, 329, 331, 333, 335)
Z_DZ = 0.03125
Z_PROBES = [-10.0 + (j + 0.5) * Z_DZ for j in Z_JS]


def theory(sigma, lam):
    omega = 2.0 * np.pi * C / lam
    k = omega * np.sqrt(MU0 * EPS0) * np.sqrt(1.0 - 1j * sigma / (omega * EPS0))
    return k, -1.0 / k.imag


def dft(t, y, w, t0, t1):
    m = (t >= t0) & (t <= t1)
    return np.sum(y[m] * np.exp(-1j * w * t[m])) * (t[1] - t[0])


def fit_delta(zs, amps):
    zs = np.asarray(zs)
    amps = np.asarray(amps)
    _, delta_th = theory(ARGS.sigma, ARGS.lam)
    mask = zs >= getattr(ARGS, 'fit_start', 1.5) * delta_th
    if mask.sum() < 4:
        mask = np.ones_like(zs, dtype=bool)
    coef, res, _rank, _sv, _rc = np.polyfit(zs[mask], np.log(amps[mask]), 1, full=True)
    delta_meas = -1.0 / coef[0]
    resid = float(np.sqrt(res[0] / max(mask.sum() - 2, 1))) if len(res) else 0.0
    return delta_meas, resid, mask


def selftest():
    """Synthetic signal with a known delta: the fit must recover it exactly."""
    _, delta = theory(ARGS.sigma, ARGS.lam)
    w = 2.0 * np.pi * C / ARGS.lam
    t = np.arange(0.0, 60e-9, 1e-12)
    amps = []
    for z in Z_PROBES:
        env = np.exp(-((t - (16e-9 + (z + 5.0) / C)) / 2e-9) ** 2)
        amps.append(abs(dft(t, env * np.cos(w * t) * np.exp(-z / delta), w, ARGS.t0, ARGS.t1)))
    d_meas, resid, mask = fit_delta(Z_PROBES, amps)
    print(f"selftest: delta_meas = {d_meas:.6f} m   theory {delta:.6f} m   "
          f"rel {d_meas / delta - 1:+.2e}   resid {resid:.2e}")
    ok = abs(d_meas / delta - 1) < 1e-6
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return 0 if ok else 1


def main():
    global ARGS
    ap = argparse.ArgumentParser()
    ap.add_argument("diags", nargs="?", help="directory with q0.txt..q6.txt")
    ap.add_argument("--sigma", type=float, default=0.1)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--t0", type=float, default=24e-9)
    ap.add_argument("--t1", type=float, default=50e-9)
    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--dz", type=float, default=0.03125, help="cell size (1280^2 uses 0.015625)")
    ap.add_argument("--j0", type=int, default=323, help="first probe cell index")
    ap.add_argument("--dj", type=int, default=2, help="probe cell index step")
    ap.add_argument("--fit-start", type=float, default=1.5,
                    help="start the ln|E| fit at this many skin depths")
    ap.add_argument("--selftest", action="store_true")
    ARGS = ap.parse_args()

    global Z_PROBES
    Z_PROBES = [-10.0 + (ARGS.j0 + k * ARGS.dj + 0.5) * ARGS.dz for k in range(7)]

    if ARGS.selftest:
        return selftest()

    k, delta_th = theory(ARGS.sigma, ARGS.lam)
    w = 2.0 * np.pi * C / ARGS.lam
    amps = []
    for i, z in enumerate(Z_PROBES):
        d = np.loadtxt(f"{ARGS.diags}/p{i}.txt", skiprows=1)
        amps.append(abs(dft(d[:, 1], d[:, 6], w, ARGS.t0, ARGS.t1)))

    d_meas, resid, mask = fit_delta(Z_PROBES, amps)
    rel = d_meas / delta_th - 1.0

    print(f"k = {k.real:.6f} {k.imag:+.6f}i   delta_theory = {delta_th:.6f} m")
    print(f"probe z / |E| (ln) :")
    for z, a in zip(Z_PROBES, amps):
        print(f"   z={z:7.4f}  |E|={a:.6e}  ln|E|={np.log(a):+.6f}"
              + ("   <- excluded (z < 1.5 delta)" if not mask[Z_PROBES.index(z)] else ""))
    print(f"slope fit ({mask.sum()} points, z >= {ARGS.fit_start} delta): delta_meas = {d_meas:.6f} m   "
          f"rel {rel:+.3%}   resid {resid:.2e}")
    print("PASS" if abs(rel) <= ARGS.tol else "FAIL")
    return 0 if abs(rel) <= ARGS.tol else 1


if __name__ == "__main__":
    raise SystemExit(main())
