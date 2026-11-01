#!/usr/bin/env python3
"""Convergence order of the PEC-cylinder scattering against the 2D Mie series.

For every (polarization, resolution) run:

  1. read the last 8 full plotfiles (one per T/8) of the main field component,
  2. extract the time-harmonic amplitude with an 8-point DFT at omega,
  3. align one global phase to the analytic solution (plane-wave reference),
  4. near field: L2 relative error over the annulus 1.5R <= r <= 2.5R,
  5. RCS: project the scattered field on the circle r = 2R onto cylindrical
     harmonics, propagate each harmonic to the far field and compare sigma(phi)
     with the Mie result (L2 relative error),
  6. fit log(err) vs log(dx) to get the observed order.

Usage:  python3 analyze_mie_convergence.py <root> [pol ...]
        <root> holds sub-directories <pol>_<n>/ with diags/ and inputs.
"""
import glob
import math
import os
import re
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import hankel1

import mie_cylinder as mie

C = 299792458.0
LAM = 1.0
K = 2.0 * math.pi / LAM
OMEGA = K * C
T = LAM / C
R = 1.0
R_IN, R_OUT = 1.5 * R, 2.5 * R
R_RCS = 2.0 * R
FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def read_pf(path):
    """Read a 2D single-box plotfile; return the component array (nz, nx)."""
    hdr = open(os.path.join(path, "Header")).read().split("\n")
    nv = int(hdr[1])
    boxline = [l for l in hdr if l.startswith("((")][0]
    nums = re.findall(r"\((\d+),(\d+)\)", boxline)
    lo = (int(nums[0][0]), int(nums[0][1]))
    hi = (int(nums[1][0]), int(nums[1][1]))
    nx, nz = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1
    raw = open(os.path.join(path, "Level_0", "Cell_D_00000"), "rb").read()
    need = nx * nz * nv * 8
    arr = np.frombuffer(raw[len(raw) - need:], dtype="<f8").reshape(nv, nz, nx)
    return arr[0]


def load_frames(run_dir, interval):
    diags = sorted(glob.glob(os.path.join(run_dir, "diags", "d1*")),
                   key=lambda p: int(re.search(r"d1(\d+)", p).group(1)))
    steps = [int(re.search(r"d1(\d+)", d).group(1)) for d in diags]
    for i in range(len(diags) - 8, -1, -1):
        sel = list(range(i, i + 8))
        if all(steps[j + 1] - steps[j] == interval for j in sel[:-1]):
            return [diags[j] for j in sel]
    return diags[-8:]


def grid_from_inputs(run_dir):
    txt = open(os.path.join(run_dir, "inputs")).read()
    n = re.search(r"^amr\.n_cell\s*=\s*(\d+)\s+(\d+)", txt, re.M)
    lo = re.search(rf"^geometry\.prob_lo\s*=\s*({FLOAT})\s+({FLOAT})", txt, re.M)
    hi = re.search(rf"^geometry\.prob_hi\s*=\s*({FLOAT})\s+({FLOAT})", txt, re.M)
    nx, nz = int(n.group(1)), int(n.group(2))
    x0, z0 = float(lo.group(1)), float(lo.group(2))
    x1, z1 = float(hi.group(1)), float(hi.group(2))
    x = x0 + (np.arange(nx) + 0.5) * (x1 - x0) / nx
    z = z0 + (np.arange(nz) + 0.5) * (z1 - z0) / nz
    return x, z


def pf_time(f):
    """Physical time stored on the 5th line of a plotfile Header."""
    with open(os.path.join(f, 'Header')) as fh:
        for _ in range(4):
            fh.readline()
        return float(fh.readline())


def harmonic_amplitude(frames, x, z):
    # The frame spacing in steps is an integer (d1.intervals), so it is NOT
    # exactly T/8.  Using the physical time of each frame (instead of j*T/8)
    # removes the resulting phase error, which was ~34 deg at 40^2 (interval=2).
    acc = np.zeros((len(z), len(x)), dtype=complex)
    for f in frames:
        acc += read_pf(f) * np.exp(1j * OMEGA * pf_time(f))
    return 2.0 / len(frames) * acc


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    extra = sys.argv[2:]
    vac_root = None
    if "--vacuum" in extra:
        i = extra.index("--vacuum")
        vac_root = extra[i + 1]
        extra = extra[:i] + extra[i + 2:]
    pols = extra or ["ey", "by"]

    summary = {}
    for pol in pols:
        runs = sorted(glob.glob(os.path.join(root, f"{pol}_*")),
                      key=lambda p: int(p.split("_")[-1]))
        results = []
        for run in runs:
            n = int(run.split("_")[-1])
            txt = open(os.path.join(run, "inputs")).read()
            interval = int(re.search(r"^d1\.intervals\s*=\s*(\d+)", txt, re.M).group(1))
            frames = load_frames(run, interval)
            if len(frames) < 8:
                print(f"{pol}_{n}: only {len(frames)} frames, skipped")
                continue

            x, z = grid_from_inputs(run)
            dx = x[1] - x[0]
            A = harmonic_amplitude(frames, x, z)
            # normalize by the incident amplitude, otherwise the L2 error is
            # dominated by E0 (=1e6 V/m) instead of the discretization error
            _txt = open(os.path.join(run, 'inputs')).read()
            for _line in _txt.splitlines():
                if _line.strip().startswith('my_constants.E0'):
                    A = A / float(_line.split('=')[1].strip())
                    break
            X, Z = np.meshgrid(x, z, indexing="xy")
            r = np.hypot(X, Z)
            th = np.arctan2(X, Z)

            Am = mie.field(K, R, r, th, pol)
            Am = np.where(r > 1e-9, Am, 0.0)      # H_n(0) is singular at the origin
            mask = (r >= R_IN) & (r <= R_OUT)

            phase = np.angle(np.sum(A[mask] * np.conj(Am[mask])))
            A_al = A * np.exp(1j * phase)

            # Near field: compare *scattered* fields when a numerical incident field is
            # available, otherwise compare total fields (which keeps the ~14% floor).
            _A_vac_tot = None
            if vac_root is not None:
                _vd = os.path.join(vac_root, "%s_%d" % (pol, n))
                _tf = os.path.join(_vd, "inputs")
                _iv = int(re.search(r"^d1\.intervals\s*=\s*(\d+)", open(_tf).read(), re.M).group(1))
                _vf = load_frames(_vd, _iv)
                _A_vac_tot = harmonic_amplitude(_vf, x, z)
                for _line in open(_tf).read().splitlines():
                    if _line.strip().startswith("my_constants.E0"):
                        _A_vac_tot = _A_vac_tot / float(_line.split("=")[1].strip())
                        break
            if _A_vac_tot is not None:
                _A_inc = _A_vac_tot * np.exp(1j * phase)
                _A_sc_num = A_al - _A_inc
                _A_sc_th = Am - np.exp(1j * K * Z)
                err_nf = math.sqrt(np.sum(np.abs(_A_sc_num[mask] - _A_sc_th[mask]) ** 2)
                                   / np.sum(np.abs(_A_sc_th[mask]) ** 2))
            else:
                err_nf = math.sqrt(np.sum(np.abs(A_al[mask] - Am[mask]) ** 2)
                                   / np.sum(np.abs(Am[mask]) ** 2))

            # RCS: project the scattered field on r = R_RCS onto harmonics
            if vac_root is not None:
                Asc = A_al - _A_vac_tot * np.exp(1j * phase)
            else:
                Asc = A_al - np.exp(1j * K * Z)
            th_c = np.linspace(0, 2 * np.pi, 721, endpoint=False)
            pts = np.stack([R_RCS * np.cos(th_c), R_RCS * np.sin(th_c)], axis=-1)
            fsc = (RegularGridInterpolator((z, x), Asc.real, bounds_error=False,
                                          fill_value=0.0)(pts)
                   + 1j * RegularGridInterpolator((z, x), Asc.imag, bounds_error=False,
                                                  fill_value=0.0)(pts))
            n_max = mie.nmax(K * R)
            ns = np.arange(-n_max, n_max + 1)
            F = np.array([np.mean(fsc * np.exp(-1j * n_ * th_c)) for n_ in ns])
            cn = -F / ((1j ** ns) * hankel1(ns, K * R_RCS))
            sig_sim = (4.0 / K) * np.abs(
                np.sum((1j ** ns)[:, None] * cn[:, None]
                       * np.exp(1j * ns[:, None] * th_c[None, :]), axis=0)) ** 2
            sig_mie = mie.rcs(K, R, th_c, pol)
            err_rcs = math.sqrt(np.sum((sig_sim - sig_mie) ** 2) / np.sum(sig_mie ** 2))

            results.append((dx, err_nf, err_rcs))
            print(f"{pol}_{n:4d}: dx={dx:.6g}  L2(near)={err_nf:.4e}  "
                  f"L2(RCS)={err_rcs:.4e}  phase={math.degrees(phase):+.1f} deg",
                  flush=True)

        if len(results) >= 2:
            dx = np.array([r[0] for r in results])
            p_nf = np.polyfit(np.log(dx), np.log([r[1] for r in results]), 1)[0]
            p_rc = np.polyfit(np.log(dx), np.log([r[2] for r in results]), 1)[0]
            print(f"  {pol.upper()}: observed order  near field = {p_nf:.3f}  "
                  f"RCS = {p_rc:.3f}")
            summary[pol] = (p_nf, p_rc)

    if summary:
        print("\nsummary (observed convergence order, log-log fit):")
        for pol, (p1, p2) in summary.items():
            print(f"  {pol.upper()}: near field {p1:.2f}, RCS {p2:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
