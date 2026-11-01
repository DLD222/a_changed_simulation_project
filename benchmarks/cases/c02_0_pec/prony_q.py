#!/usr/bin/env python3
"""Multi-mode complex-exponential fit (matrix-pencil / Prony) for resonator Q.

The envelope-slope method fails whenever the probe sees a standing wave (several
modes plus reflections): |Ey(t)| oscillates and ln|Ey| is not linear. This script
instead decomposes the *complex* time series into damped exponentials

    x(t) ~ sum_k A_k exp((-gamma_k + i 2 pi f_k) t)

so every mode comes out at once with its own frequency and decay. For each mode

    Q_k = pi f_k / gamma_k

Method: build the analytic signal with a Hilbert transform (the recorded field is
real), form two Hankel matrices from the samples and solve the generalised
eigenvalue problem (matrix pencil). No assumption of a single mode.

Usage:
  python3 prony_q.py --probe /tmp/t7q2/diags/reducedfiles/p.txt \
      --t0 3.6e-6 --t1 6.8e-6 --order 24
"""
import argparse

import numpy as np
from scipy.signal import hilbert


def matrix_pencil(x, dt, order):
    """Return (freq, decay, amplitude, pole) for the strongest `order` modes."""
    n = len(x)
    if 2 * order > n:
        raise SystemExit(f"need 2*order <= len(x) (order={order}, n={n})")
    # Hankel matrices Y0 = x[i+j], Y1 = x[i+j+1]
    cols = n - 2 * order + 1
    Y0 = np.empty((order, cols), dtype=complex)
    Y1 = np.empty((order, cols), dtype=complex)
    for j in range(cols):
        Y0[:, j] = x[j:j + order]
        Y1[:, j] = x[j + 1:j + order + 1]

    # least-squares Prony: poles z_k = exp((-gamma + i w) dt) are the eigenvalues of
    # Y1 @ pinv(Y0)  (a square order x order problem; eig(Y1, Y0) would need Y0 square)
    poles = np.linalg.eigvals(Y1 @ np.linalg.pinv(Y0))
    poles = poles[np.isfinite(poles) & (np.abs(poles) > 0)]

    # Residues from a global Vandermonde are hopelessly ill-conditioned here
    # (n ~ 1e5, |z| ~ 1), so instead rank the modes by how close |z| is to 1:
    # the slowest-decaying pole is the resonance of interest.
    freq = np.angle(poles) / (2 * np.pi * dt)
    decay = -np.log(np.abs(poles)) / dt
    keep = (freq > 0) & (decay > 0) & (np.abs(poles) < 1.0)
    freq, decay, poles = freq[keep], decay[keep], poles[keep]
    order_by_slow = np.argsort(decay)
    return freq[order_by_slow], decay[order_by_slow], np.abs(poles[order_by_slow]), poles[order_by_slow]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", required=True)
    ap.add_argument("--col", type=int, default=6, help="column index of Ey (default 6)")
    ap.add_argument("--t0", type=float, required=True)
    ap.add_argument("--t1", type=float, required=True)
    ap.add_argument("--order", type=int, default=24)
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--dec", type=int, default=20, help="decimation factor (conditioning)")
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--c", type=float, default=299792458.0)
    args = ap.parse_args()

    d = np.loadtxt(args.probe, skiprows=1)
    t = d[:, 1]
    sig = d[:, args.col]
    m = (t >= args.t0) & (t <= args.t1)
    t, sig = t[m], sig[m]
    t, sig = t[::args.dec], sig[::args.dec]
    dt = t[1] - t[0]
    print(f"{len(sig)} samples, dt={dt:.4e} s, span {(t[-1]-t[0])*1e6:.3f} us")
    print(f"|x| max={np.abs(sig).max():.3e}, mean={sig.mean():.3e}  (mean is removed)")

    x = hilbert(sig - sig.mean())
    freq, decay, amps, _ = matrix_pencil(x, dt, args.order)

    print(f"\n{'f (Hz)':>14} {'f_meep':>9} {'gamma (1/s)':>13} {'Q':>9} {'|z|':>9}")
    for f, g, A in zip(freq[:args.top], decay[:args.top], amps[:args.top]):
        fm = f * args.a / args.c
        print(f"{f:14.5e} {fm:9.5f} {g:13.4e} {np.pi*f/g:9.1f} {abs(A):9.5f}")


if __name__ == "__main__":
    main()
