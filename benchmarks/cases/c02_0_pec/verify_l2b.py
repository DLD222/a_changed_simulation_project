#!/usr/bin/env python3
"""L2b: effective wall position of a grid-aligned PEC plane, ECT vs Yee.

The plane sits exactly on the mesh (z = 0). Such a surface still creates "cut"
cells in AMReX: the open cell that touches the plane has volume fraction 1 but
is flagged as non-regular, and the staircase masking used by the Yee solver
treats it as closed. The conformal geometry used by ECT keeps it open, so the
two solvers place the reflecting wall half a cell apart.

Criterion (revised after the M2A code audit, see M2A_CODE_AUDIT.md P1):

  * the **round-trip time** ``t_ref - t_inc`` is used, not the absolute arrival
    time. The probe reads a carrier peak, whose phase drifts with ``dz`` and with
    the sampling, so ``t_ref`` alone is not a robust quantity at higher
    resolution - both ``t_inc`` and ``t_ref`` move together and cancel;
  * ``dz`` is parsed from the run's own ``warpx_used_inputs``/``inputs``
    (``amr.n_cell`` + ``geometry.prob_lo/hi``); the old hard-coded 0.25 default
    silently gave the wrong analytic value for any other resolution;
  * the peak ratio ``r`` is reported for information only. It is a carrier-peak
    ratio, not a proper reflection coefficient, so it is not used as a gate.

Acceptance:
  * |round_trip(ECT) - 2*(z_probe - z_wall)/c| <= tol (default 0.5 ns),
  * ECT is closer to that analytic round-trip time than Yee.

Usage:  python3 verify_l2b.py <ect_dir> <yee_dir> [dz] [probe_z] [pulse_z] [z_wall] [tol_ns]
Any argument left as "-" is derived from the run inputs / probe file.
"""
import os
import re
import sys

C = 299792458.0
FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("[") or s.startswith("#"):
                continue
            rows.append([float(v) for v in line.split()])
    return rows


def parse_inputs(run_dir):
    """Return (dz, z_pulse) parsed from the run's inputs, or (None, None)."""
    text = None
    for name in ("warpx_used_inputs", "inputs"):
        path = os.path.join(run_dir, name)
        if os.path.exists(path):
            text = open(path).read()
            break
    if text is None:
        return None, None

    dz = None
    n = re.search(r"^amr\.n_cell\s*=\s*(\d+)\s+(\d+)", text, re.M)
    lo = re.search(rf"^geometry\.prob_lo\s*=\s*({FLOAT})\s+({FLOAT})", text, re.M)
    hi = re.search(rf"^geometry\.prob_hi\s*=\s*({FLOAT})\s+({FLOAT})", text, re.M)
    if n and lo and hi:
        nz = int(n.group(2))
        dz = (float(hi.group(2)) - float(lo.group(2))) / nz

    zc = re.search(rf"^my_constants\.zc\s*=\s*({FLOAT})", text, re.M)
    z_pulse = float(zc.group(1)) if zc else None
    return dz, z_pulse


def peaks(path):
    """Return (t_inc, Ey_inc, t_ref, Ey_ref) using the analytic-free heuristic:
    incident = largest |Ey| in the first half, reflected = largest |Ey| after it.

    The caller anchors the windows with the analytic times, so the heuristic is
    only used to locate the two pulses inside their windows.
    """
    rows = load(path)
    return rows


def peak_in(rows, t_lo, t_hi):
    best = None
    for r in rows:
        if t_lo <= r[1] <= t_hi and (best is None or abs(r[6]) > abs(best[1])):
            best = (r[1], r[6])
    return best


def measure(run_dir, dz, z_pulse, z_wall):
    path = os.path.join(run_dir, "diags", "reducedfiles", "fp1.txt")
    rows = load(path)
    z_probe = rows[0][4]

    # Analytic anchors, wide enough to contain either wall convention
    t_inc_exact = (z_probe - z_pulse) / C
    t_ref_wall = (z_wall - z_pulse + z_wall - z_probe) / C
    span = t_ref_wall - t_inc_exact          # analytic round trip (positive)
    half = 0.4 * span

    inc = peak_in(rows, t_inc_exact - half, t_inc_exact + half)
    ref = peak_in(rows, t_ref_wall - half, t_ref_wall + half)
    if inc is None or ref is None:
        return None
    round_trip = ref[0] - inc[0]
    return {
        "z_probe": z_probe,
        "t_inc": inc[0],
        "ey_inc": inc[1],
        "t_ref": ref[0],
        "ey_ref": ref[1],
        "round_trip": round_trip,
        "r": ref[1] / inc[1],
    }


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    ect_dir, yee_dir = sys.argv[1], sys.argv[2]
    dz_arg = sys.argv[3] if len(sys.argv) > 3 else "-"
    probe_arg = sys.argv[4] if len(sys.argv) > 4 else "-"
    pulse_arg = sys.argv[5] if len(sys.argv) > 5 else "-"
    wall_arg = sys.argv[6] if len(sys.argv) > 6 else "-"
    tol_ns = float(sys.argv[7]) if len(sys.argv) > 7 else 0.5

    dz_parsed, pulse_parsed = parse_inputs(ect_dir)
    dz = float(dz_arg) if dz_arg != "-" else (dz_parsed if dz_parsed else 0.25)
    z_pulse = float(pulse_arg) if pulse_arg != "-" else (pulse_parsed if pulse_parsed is not None else -5.0)
    z_wall = float(wall_arg) if wall_arg != "-" else 0.0

    ect = measure(ect_dir, dz, z_pulse, z_wall)
    yee = measure(yee_dir, dz, z_pulse, z_wall)
    if ect is None or yee is None:
        print("could not locate both pulses; check the probe windows")
        return 1

    z_probe = ect["z_probe"] if probe_arg == "-" else float(probe_arg)
    t_inc_exact = (z_probe - z_pulse) / C
    rt_exact = (z_wall - z_pulse + z_wall - z_probe) / C - t_inc_exact
    rt_inside = rt_exact - 2.0 * dz / C

    print(f"dz = {dz:.6g} m (parsed from inputs: {'yes' if dz_parsed else 'no, fallback'}), "
          f"z_probe = {z_probe:+.4f} m, z_pulse = {z_pulse:+.4f} m, z_wall = {z_wall:+.4f} m")
    print(f"analytic round trip: {rt_exact * 1e9:.3f} ns (wall at z_wall), "
          f"{rt_inside * 1e9:.3f} ns (wall one cell inside)")

    for tag, m in (("ECT", ect), ("YEE", yee)):
        print(f"{tag}: t_inc={m['t_inc'] * 1e9:8.3f} ns (Ey={m['ey_inc']:+.4e})  "
              f"t_ref={m['t_ref'] * 1e9:8.3f} ns (Ey={m['ey_ref']:+.4e})  "
              f"round_trip={m['round_trip'] * 1e9:7.3f} ns  "
              f"err={m['round_trip'] * 1e9 - rt_exact * 1e9:+6.3f} ns  r={m['r']:+.4f}")

    err_ect = abs(ect["round_trip"] - rt_exact)
    err_yee = abs(yee["round_trip"] - rt_exact)

    checks = [
        (f"|round_trip(ECT) - analytic| = {err_ect * 1e9:.3f} ns <= {tol_ns} ns",
         err_ect <= tol_ns * 1e-9),
        (f"ECT closer to the analytic round trip than Yee "
         f"({err_ect * 1e9:.3f} vs {err_yee * 1e9:.3f} ns)", err_ect < err_yee),
    ]
    ok = True
    for label, good in checks:
        ok = ok and good
        print(f"  {'OK  ' if good else 'FAIL'} {label}")

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
