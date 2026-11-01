#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Oblique-incidence Fresnel R/T by the Poynting-flux method (no k needed).

Why not two-probe wave decomposition: on the incident side the field is
incident+reflected, so `measure_k`'s single-forward-wave assumption fails and k1
comes out ~12% wrong, which wrecks the 2x2 solve. Because the wave is uniform in
x (periodic, fixed kx), the z-component of the Poynting flux on a z=const plane is
enough:

    S_z = (Ex*By - Ey*Bx)/mu0
    R = (W_vac - W_med)/W_vac        (incident-side probe)
    T = W_med / W_vac                (transmission-side probe)

The vacuum reference (same grid/source, no interface) makes both runs share the
same Yee dispersion, so no analytic k or omega enters.

Polarisations (k_hat = (sin t, 0, cos t)):
  s: E = E0 yhat          -> Ey,  Bx = -cos t * Ey/c,  Bz = +sin t * Ey/c
  p: E = E0(cos t,0,-sin t) -> Ex, Ez,  By = +E0/c

Usage:
  python3 verify_oblique_rt.py /tmp/t2p --theta 30 --eps2 4 --pol p --run
"""
import argparse
import math
import os
import subprocess

import numpy as np

C = 299792458.0
MU0 = 1.25663706212e-6
BIN = "build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB"


def make_inputs(root, args, pol, interfaced):
    """Write <root>/<med|vac>/inputs."""
    d = os.path.join(root, "med" if interfaced else "vac")
    os.makedirs(d, exist_ok=True)

    n, L = args.ncell, args.L
    k0 = 2 * math.pi / args.lam
    th = math.radians(args.theta)
    kx, kz = k0 * math.sin(th), k0 * math.cos(th)
    w, z0 = args.w, args.z0

    def ph(comp):
        return f"exp(-((z-({z0}))/{w})^2)*cos({kx:.10f}*x+{kz:.10f}*z)"

    if pol == "s":
        ex, ey, ez = "0", ph("y"), "0"
        bx = f"-{math.cos(th):.10f}*" + ph("y") + "/299792458.0"
        by, bz = "0", f"{math.sin(th):.10f}*" + ph("y") + "/299792458.0"
    else:
        ex = f"{math.cos(th):.10f}*" + ph("x")
        ey = "0"
        ez = f"-{math.sin(th):.10f}*" + ph("z")
        bx, bz = "0", "0"
        by = ph("y") + "/299792458.0"

    eps2 = args.eps2 if interfaced else 1.0
    lines = [
        f"max_step = {args.nsteps}",
        f"amr.n_cell = {n} {n}",
        f"amr.max_grid_size = {n}",
        "amr.max_level = 0",
        "geometry.dims = 2",
        f"geometry.prob_lo = {-L/2} {-L/2}",
        f"geometry.prob_hi =  {L/2}  {L/2}",
        "warpx.cfl = 0.4",
        "algo.maxwell_solver = yee",
        "algo.particle_shape = 1",
        "algo.em_solver_medium = macroscopic",
        "algo.macroscopic_subpixel = 1",
        "macroscopic.epsilon = 1.0*epsilon0",
        f"macroscopic.epsilon2 = {eps2}*epsilon0",
        "macroscopic.mu = mu0",
        "macroscopic.sigma = 0",
    ]
    if interfaced:
        lines.append('macroscopic.interface_function(x,y,z) = "-z"')
    lines += [
        "boundary.field_lo = periodic pec",
        "boundary.field_hi = periodic pec",
        'warpx.eb_implicit_function = "z - 100"',
        "warpx.E_ext_grid_init_style = parse_E_ext_grid_function",
        f'warpx.Ex_external_grid_function(x,y,z) = "{ex}"',
        f'warpx.Ey_external_grid_function(x,y,z) = "{ey}"',
        f'warpx.Ez_external_grid_function(x,y,z) = "{ez}"',
        "warpx.B_ext_grid_init_style = parse_B_ext_grid_function",
        f'warpx.Bx_external_grid_function(x,y,z) = "{bx}"',
        f'warpx.By_external_grid_function(x,y,z) = "{by}"',
        f'warpx.Bz_external_grid_function(x,y,z) = "{bz}"',
    ]
    lines.append("warpx.reduced_diags_names = q1 q3")
    for tag, z in (("q1", args.z_inc), ("q3", args.z_trans)):
        lines += [f"{tag}.type = FieldProbe",
                  f"{tag}.intervals = 1",
                  f"{tag}.probe_geometry = Point",
                  f"{tag}.x_probe = 0.0",
                  f"{tag}.z_probe = {z}"]
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(d, "inputs"), "w") as f:
        f.write(txt)
    return d


def flux(path, tag, t1, t2):
    x = np.loadtxt(os.path.join(path, "diags/reducedfiles", tag + ".txt"), skiprows=1)
    t, ex, ey, bx, by = x[:, 1], x[:, 5], x[:, 6], x[:, 8], x[:, 9]
    m = (t >= t1) & (t <= t2)
    return float(np.sum((ex[m] * by[m] - ey[m] * bx[m]) / MU0) * (t[1] - t[0]))


def fresnel(theta_deg, eps2, pol):
    n1, n2 = 1.0, math.sqrt(eps2)
    t1 = math.radians(theta_deg)
    s2 = n1 * math.sin(t1) / n2
    if s2 >= 1.0:
        return 1.0, 0.0
    t2 = math.asin(s2)
    c1, c2 = math.cos(t1), math.cos(t2)
    if pol == "s":
        r = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
        t = 2 * n1 * c1 / (n1 * c1 + n2 * c2)
    else:
        r = (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2)
        t = 2 * n1 * c1 / (n2 * c1 + n1 * c2)
    T = (n2 * c2) / (n1 * c1) * t * t
    return r * r, T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--theta", type=float, default=30.0)
    ap.add_argument("--eps2", type=float, default=4.0)
    ap.add_argument("--pol", choices=("s", "p"), default="s")
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--ncell", type=int, default=480)
    ap.add_argument("--L", type=float, default=24.0)
    ap.add_argument("--w", type=float, default=2.5, help="envelope half-width (>=3 lambda_z)")
    ap.add_argument("--z0", type=float, default=-5.0)
    ap.add_argument("--nsteps", type=int, default=1300)
    ap.add_argument("--z-inc", type=float, default=-2.0)
    ap.add_argument("--z-trans", type=float, default=0.4)
    ap.add_argument("--tmax", type=float, default=58e-9)
    ap.add_argument("--make", action="store_true")
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    if args.make:
        for iff in (False, True):
            make_inputs(args.root, args, args.pol, iff)
        print("wrote", args.root)

    if args.run:
        for sub in ("vac", "med"):
            d = os.path.join(args.root, sub)
            with open(os.path.join(d, "run.log"), "w") as f:
                r = subprocess.run([BIN, "inputs"], cwd=d, stdout=f, stderr=subprocess.STDOUT,
                                   env={**os.environ, "OMP_NUM_THREADS": "8"})
            print(f"  ran {sub}: exit={r.returncode}")

    vac, med = os.path.join(args.root, "vac"), os.path.join(args.root, "med")
    Wi = flux(vac, "q1", 0.0, args.tmax)
    Wr = flux(med, "q1", 0.0, args.tmax)
    R = (Wi - Wr) / Wi
    T = flux(med, "q3", 0.0, args.tmax) / flux(vac, "q3", 0.0, args.tmax)
    Rth, Tth = fresnel(args.theta, args.eps2, args.pol)

    print(f"pol={args.pol} theta={args.theta} eps2={args.eps2} "
          f"ncell={args.ncell} (ppw~{args.lam/(args.L/args.ncell):.0f})")
    print(f"  R = {R:.4f}   theory {Rth:.4f}   rel {R/Rth-1:+.2%}")
    print(f"  T = {T:.4f}   theory {Tth:.4f}   rel {T/Tth-1:+.2%}")
    print(f"  R+T = {R+T:.4f}")
    if args.pol == "p":
        tb = math.degrees(math.atan(math.sqrt(args.eps2)))
        print(f"  Brewster angle (this eps2) = {tb:.2f} deg")


if __name__ == "__main__":
    main()
