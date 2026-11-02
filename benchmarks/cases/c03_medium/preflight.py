#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Pre-flight checks for a macroscopic (c03) run.

Motivation: two real mistakes cost us a full debugging cycle each
  * `macroscopic.sigma = 0.01` (global scalar) made the "vacuum" side lossy, so the
    incident wave was attenuated by e^{-43} BEFORE reaching the interface and the
    measured transmission was ~0 (see MEDIUM_MULTIGRID_ROOTCAUSE_AND_FIX.md 5.0x).
  * a fixed-index dump read a FAB whose valid box did not contain that index and
    returned shifted garbage (same section).

Run this BEFORE interpreting any R / T number:

    python3 preflight.py /tmp/case/inputs
    python3 preflight.py /tmp/case/inputs --zc -4.0 --probe -2.0 1.0

Exit code is 1 when a check fails, so it can gate a script.
"""
import argparse
import math
import re
import sys

C = 299792458.0
EPS0 = 8.8541878188e-12
MU0 = 1.25663706212e-06


def read_inputs(path):
    kv = {}
    for raw in open(path, encoding="utf-8", errors="replace"):
        line = raw.split("#", 1)[0].strip()
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        kv[k.strip()] = v.strip().strip('"')
    return kv


def as_float(expr, consts):
    """Evaluate a simple numeric / my_constants expression."""
    e = expr.strip().strip('"')
    try:
        return float(eval(e, {"__builtins__": {}}, dict(consts)))
    except Exception:
        return None


def sigma_profile(kv, consts):
    """Return (kind, sigma1, sigma2) where sigma2 applies to z>0 if a step is used."""
    if "macroscopic.sigma_function(x,y,z)" in kv:
        expr = kv["macroscopic.sigma_function(x,y,z)"]
        m = re.match(r"\s*([\w.*+\-/ ]+?)\s*\*\s*\(?\s*z\s*>\s*([-\d.eE+]+)\s*\)?", expr)
        if m:
            amp = as_float(m.group(1), consts)
            return ("step", 0.0 if amp is None else 0.0, amp, float(m.group(2)))
        return ("function", None, None, 0.0)
    if "macroscopic.sigma" in kv:
        s = as_float(kv["macroscopic.sigma"], consts)
        return ("global", s, s, 0.0)
    return ("none", 0.0, 0.0, 0.0)


def eps_profile(kv, consts):
    """Return (eps1, eps2, eps_yy2) for s-polarisation."""
    eps1 = as_float(kv.get("macroscopic.epsilon", "eps0"), consts) or EPS0
    eps2 = eps1
    if "macroscopic.epsilon2" in kv:
        v = as_float(kv["macroscopic.epsilon2"], consts)
        if v:
            eps2 = v
    # diagonal tensor input overrides the second material
    k = "macroscopic.epsilon2_tensor_yy"
    if k in kv:
        v = as_float(kv[k], consts)
        if v:
            eps2 = v
    return eps1, eps2


def complex_k(omega, eps, sigma):
    return omega * math.sqrt(MU0 * eps) * cmath_sqrt(1.0 - 1j * sigma / (omega * eps))


def cmath_sqrt(z):
    r = math.hypot(z.real, z.imag)
    re = math.sqrt(max(0.0, 0.5 * (r + z.real)))
    im = math.sqrt(max(0.0, 0.5 * (r - z.real)))
    if z.imag < 0:
        im = -im
    return complex(re, im)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs")
    ap.add_argument("--zc", type=float, default=None, help="source z (default: my_constants.zc)")
    ap.add_argument("--probe", type=float, nargs="+", default=[-2.0, 0.05, 1.0])
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    kv = read_inputs(a.inputs)
    consts = {"eps0": EPS0, "mu0": MU0, "c": C, "pi": math.pi}
    for k, v in kv.items():
        if k.startswith("my_constants."):
            f = as_float(v, consts)
            if f is not None:
                consts[k.split(".", 1)[1]] = f

    lam = as_float(kv.get("my_constants.lam", "1.0"), consts) or 1.0
    omega = 2.0 * math.pi * C / lam
    lo = [float(x) for x in kv.get("geometry.prob_lo", "-1 -1").split()[:2]]
    hi = [float(x) for x in kv.get("geometry.prob_hi", "1 1").split()[:2]]
    ncell = [int(x) for x in kv.get("amr.n_cell", "64 64").split()[:2]]
    mgs = int(as_float(kv.get("amr.max_grid_size", "64"), consts) or 64)
    dz = (hi[1] - lo[1]) / ncell[1]
    dx = (hi[0] - lo[0]) / ncell[0]

    zc = a.zc if a.zc is not None else (as_float(kv.get("my_constants.zc", "nan"), consts))
    has_iface = "macroscopic.interface_function(x,y,z)" in kv

    eps1, eps2 = eps_profile(kv, consts)
    kind, s1, s2, zstep = sigma_profile(kv, consts)

    fails = []
    print("=" * 78)
    print(f"materials   eps1 = {eps1:.6e}   eps2 = {eps2:.6e}  (eps2/eps1 = {eps2/eps1:.3f})")
    print(f"            interface: {'YES' if has_iface else 'no'}")
    print(f"            sigma: {kind}"
          + (f"  sigma1 = {s1}  sigma2 = {s2}  (step at z={zstep})" if kind == "step" else
             (f"  sigma = {s1}" if kind == "global" else "")))
    if kind == "global" and s1 not in (None, 0.0) and has_iface:
        fails.append("sigma is a GLOBAL constant but an interface is present: the 'vacuum' "
                     "side is lossy too (see 5.0x). Use sigma_function(x,y,z) = \"s*(z>0)\".")

    # ---- attenuation budget from the source to the interface -------------------
    if zc is not None and not math.isnan(zc):
        seg_len = abs(zc - zstep if has_iface else zc)
        k1 = complex_k(omega, eps1, s1 if kind != "step" else 0.0)
        d1 = 1.0 / abs(k1.imag) if abs(k1.imag) > 0 else float("inf")
        att = seg_len / d1
        print(f"path        source z={zc:.3f} -> interface z={zstep:.3f}: {seg_len:.3f} m in "
              f"material 1 (delta1 = {d1:.4f} m)")
        print(f"            e^(-int dz/delta1) = {math.exp(-att):.3e}"
              + ("   <-- incident never reaches the interface!" if att > 5 else ""))
        if att > 5:
            fails.append(f"incident is attenuated by e^-{att:.1f} before the interface "
                         f"(material 1 has sigma={s1}); the R/T numbers are meaningless.")

    # ---- expected values -------------------------------------------------------
    k2 = complex_k(omega, eps2, s2 if kind == "step" else (0.0 if kind == "none" else s1))
    d2 = 1.0 / abs(k2.imag) if abs(k2.imag) > 0 else float("inf")
    eta1 = math.sqrt(MU0 / eps1)
    eta2 = complex(MU0, 0) / complex(0, 1)  # placeholder to keep types simple
    eta2 = complex(math.sqrt(MU0 / eps2), 0) / cmath_sqrt(1.0 - 1j * (s2 if kind == "step" else s1) / (omega * eps2))
    r = (eta2 - eta1) / (eta2 + eta1)
    print(f"expected    delta2 = {d2:.4f} m   |eta2| = {abs(eta2):.2f} ohm   "
          f"R = {abs(r)**2:.4f}   T(interface) = {1-abs(r)**2:.4f}")
    print(f"            k2 = {k2.real:.4f} {k2.imag:+.4f}i   lambda2 = {2*math.pi/abs(k2.real):.4f} m   "
          f"delta2/dz = {d2/dz:.2f}")

    # ---- timing ---------------------------------------------------------------
    v1 = C / math.sqrt(eps1 / EPS0)
    v2 = C / math.sqrt(eps2 / EPS0)
    if zc is not None and not math.isnan(zc):
        for zp in a.probe:
            if has_iface and zp > zstep:
                t_inc = abs(zc - zstep) / v1 + abs(zstep - zp) / v2
            else:
                t_inc = abs(zc - zp) / v1
            t_refl = abs(zc - zstep) / v1 + abs(zstep - zp) / (v1 if zp < zstep else v2)
            t_pec = (abs(zc - zstep) / v1 + abs(zstep - hi[1]) / v2 + abs(hi[1] - zp) / v2)
            print(f"timing      probe z={zp:+.4f}: incident {t_inc*1e9:6.2f} ns   "
                  f"interface echo {t_refl*1e9:6.2f} ns   PEC(+{hi[1]:.1f}) echo {t_pec*1e9:6.2f} ns")

    # ---- index / FAB validity -------------------------------------------------
    print(f"grid        dz = {dz:.6f}  dx = {dx:.6f}  max_grid_size = {mgs}  "
          f"boxes along z: {hi[1]-lo[1]:.3f}/{mgs} = {(hi[1]-lo[1])/mgs:.2f}")
    boxes = []
    nb = int(math.ceil(ncell[1] / mgs))
    for b in range(nb):
        j0 = b * mgs
        j1 = min(ncell[1] - 1, j0 + mgs - 1)
        boxes.append((j0, j1, lo[1] + j0 * dz, lo[1] + j1 * dz))
    print("            z-boxes: " + ", ".join(f"[{a0},{a1}]" for a0, a1, _, _ in boxes))
    for zp in a.probe:
        j = (zp - lo[1]) / dz
        node = abs(j - round(j)) < 1e-9
        inb = [f"[{a0},{a1}]" for a0, a1, _, _ in boxes if a0 <= round(j) <= a1]
        print(f"            probe z={zp:+.4f} -> j={j:.4f} "
              f"({'NODE' if node else 'cell-centre/interpolated'}), "
              f"in {inb[0] if inb else 'OUTSIDE'}")

    print("=" * 78)
    if fails:
        print("PRE-FLIGHT: FAIL")
        for f in fails:
            print("  * " + f)
        return 1
    print("PRE-FLIGHT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
