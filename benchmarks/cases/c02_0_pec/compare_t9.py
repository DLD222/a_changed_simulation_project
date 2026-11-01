#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""T9: compare the WarpX dielectric-cylinder scattered field with the analytic Mie series.

The WarpX run is CW (plane wave along +z), so the scattered field is obtained by
subtracting a vacuum reference run with the same source and grid:

    E_sc(x, z) = E_tot(x, z) - E_vac(x, z)

`mie_cylinder_dielectric.scattered` gives the exact near field for unit incident
amplitude, so the comparison is done on the annulus 1.5 R < r < 2.5 R after
normalising by the source amplitude E0.

Usage:
  python3 compare_t9.py --tot /tmp/t9 --vac /tmp/t9/vac --lam 1.0 --R 1.0 --m 2.0
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_mie_convergence import read_pf          # noqa: E402
import mie_cylinder_dielectric as mie                 # noqa: E402


def last_frame(d):
    ds = sorted(glob.glob(os.path.join(d, "diags", "d1*")),
                key=lambda x: int(re.search(r"d1(\d+)", x).group(1)))
    if not ds:
        raise SystemExit(f"no frames in {d}")
    return ds[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tot", required=True)
    ap.add_argument("--vac", required=True)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--R", type=float, default=1.0)
    ap.add_argument("--m", type=float, default=2.0, help="sqrt(eps2/eps1)")
    ap.add_argument("--E0", type=float, default=1.0e6)
    ap.add_argument("--rmin", type=float, default=1.5)
    ap.add_argument("--rmax", type=float, default=2.5)
    args = ap.parse_args()

    A = read_pf(last_frame(args.tot))
    B = read_pf(last_frame(args.vac))
    n = A.shape[0]
    L = 8.0                      # domain length assumed by the case template
    h = L / n
    x = -L / 2 + (np.arange(n) + 0.5) * h
    X, Z = np.meshgrid(x, x, indexing="xy")
    r = np.hypot(X, Z)
    theta = np.arctan2(X, Z)     # from +z toward +x

    S_num = (A - B) / args.E0    # unit incident amplitude
    k = 2 * np.pi / args.lam
    x_ = k * args.R
    S_th = mie.scattered(x_, args.m, k * r, theta, "ey")

    mask = (r > args.rmin * args.R) & (r < args.rmax * args.R)
    num = S_num[mask]
    th = S_th[mask]
    l2 = np.linalg.norm(num - th) / np.linalg.norm(th)
    print(f"grid {n}x{n}, dx={h:.5f}, lam={args.lam}, R={args.R}, m={args.m}")
    print(f"annulus {args.rmin}R..{args.rmax}R: {mask.sum()} points")
    print(f"  max|E_sc(num)| = {np.abs(num).max():.4e}")
    print(f"  max|E_sc(th) | = {np.abs(th).max():.4e}")
    print(f"  L2 relative error = {l2:.4f} ({l2*100:.2f}%)")
    # best-fit complex scale, to separate amplitude/phase errors from shape errors
    c = np.vdot(th, num) / np.vdot(th, th)
    print(f"  best-fit scale = {abs(c):.4f} arg={np.degrees(np.angle(c)):+.1f} deg")
    print(f"  L2 after scale = {np.linalg.norm(num - c*th)/np.linalg.norm(c*th):.4f}")


if __name__ == "__main__":
    main()
