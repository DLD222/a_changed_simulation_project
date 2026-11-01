#!/usr/bin/env python3
"""T3: Snell's law from a 2D spatial FFT of the simulated field.

For each side of the interface the transmitted field is a single plane wave, so a
2D FFT of Ey(x, z) has one dominant peak; its wavevector gives the propagation angle.
The interface normal is z, so theta = atan2(k_x, k_z).

Usage:
  python3 analyze_snell.py <case_dir> --frame-in 300 --frame-out 600 \
      --z-interface 0.0 --n1 1.0 --n2 2.0 --theta1-deg 30
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_mie_convergence import read_pf  # noqa: E402


def frames(root):
    ds = sorted(glob.glob(os.path.join(root, "diags", "d1*")),
                key=lambda x: int(re.search(r"d1(\d+)", x).group(1)))
    return [(int(re.search(r"d1(\d+)", d).group(1)), read_pf(d)) for d in ds]


def k_peak(E, dx, dz, x0, z0, mask):
    """Dominant wavevector (kx, kz) of the masked region, in rad/m."""
    nz, nx = E.shape
    x = x0 + (np.arange(nx) + 0.5) * dx
    z = z0 + (np.arange(nz) + 0.5) * dz
    X, Z = np.meshgrid(x, z, indexing="xy")
    m = mask(X, Z)
    if m.sum() < 100:
        return None
    sub = E * m
    F = np.fft.fftshift(np.fft.fft2(sub))
    kx = np.fft.fftshift(np.fft.fftfreq(nx, dx)) * 2 * np.pi
    kz = np.fft.fftshift(np.fft.fftfreq(nz, dz)) * 2 * np.pi
    i = np.unravel_index(np.argmax(np.abs(F)), F.shape)
    return kx[i[1]], kz[i[0]], np.abs(F[i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--frame-in", type=int, default=300)
    ap.add_argument("--frame-out", type=int, default=600)
    ap.add_argument("--z-interface", type=float, default=0.0)
    ap.add_argument("--n1", type=float, default=1.0)
    ap.add_argument("--n2", type=float, default=2.0)
    ap.add_argument("--theta1-deg", type=float, default=30.0)
    ap.add_argument("--x0", type=float, default=-4.0)
    ap.add_argument("--z0", type=float, default=-4.0)
    ap.add_argument("--L", type=float, default=8.0)
    args = ap.parse_args()

    fr = dict(frames(args.case))
    if args.frame_in not in fr or args.frame_out not in fr:
        print("available frames:", sorted(fr))
        raise SystemExit(1)

    E_in = fr[args.frame_in]
    nz, nx = E_in.shape
    dx = args.L / nx
    dz = args.L / nz

    # incident side: the wave has not reached the interface yet
    r1 = k_peak(E_in, dx, dz, args.x0, args.z0,
                lambda X, Z: (Z < args.z_interface - 0.3))
    # transmitted side: only after the wave crossed
    E_out = fr[args.frame_out]
    r2 = k_peak(E_out, dx, dz, args.x0, args.z0,
                lambda X, Z: (Z > args.z_interface + 0.3))

    if r1 is None or r2 is None:
        print("FFT region empty")
        raise SystemExit(1)

    kx1, kz1, a1 = r1
    kx2, kz2, a2 = r2
    th1 = np.degrees(np.arctan2(kx1, abs(kz1)))
    th2 = np.degrees(np.arctan2(kx2, abs(kz2)))

    th2_theory = np.degrees(np.arcsin(args.n1 * np.sin(np.radians(th1)) / args.n2))
    print(f"incident : kx={kx1:.4f} kz={kz1:.4f}  theta1={th1:.3f} deg  (input {args.theta1_deg:.1f})")
    print(f"transmit : kx={kx2:.4f} kz={kz2:.4f}  theta2={th2:.3f} deg")
    print(f"Snell    : n1 sin(t1) = {args.n1*np.sin(np.radians(th1)):.5f}, "
          f"n2 sin(t2) = {args.n2*np.sin(np.radians(th2)):.5f}")
    print(f"theory   : theta2 = {th2_theory:.3f} deg")
    print(f"error    : {th2 - th2_theory:+.3f} deg   (tolerance 0.5 deg)")
    print("RESULT:", "PASS" if abs(th2 - th2_theory) <= 0.5 else "FAIL")


if __name__ == "__main__":
    main()
