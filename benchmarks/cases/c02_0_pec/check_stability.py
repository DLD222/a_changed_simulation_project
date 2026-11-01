#!/usr/bin/env python3
"""Stability smoke test for the conformal (BCK-only) TE path.

Reads the FieldProbe and FieldEnergy reduced diagnostics of one run and checks

  * every probe value is finite,
  * the total field energy stays bounded: max(total)/total(t=0) <= growth_max,
  * the peak |Ey| in the second half of the run is not larger than in the first
    half (no exponential growth),
  * the energy drift at the end is within drift_max.

This is a *smoke test*: it does not probe the small-cell stability limit (that
needs a volume-fraction scan) and it cannot prove accuracy.

Usage:  python3 check_stability.py <run_dir> [growth_max] [drift_max]
"""
import os
import sys


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("[") or s.startswith("#"):
                continue
            rows.append([float(v) for v in line.split()])
    return rows


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    run_dir = sys.argv[1]
    growth_max = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
    drift_max = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05

    base = os.path.join(run_dir, "diags", "reducedfiles")
    fp = load(os.path.join(base, "fp1.txt"))
    fe = load(os.path.join(base, "fe.txt"))

    ok = True

    finite = all(v == v and abs(v) < 1e100 for r in fp for v in r)
    ok = ok and finite
    print(f"finite probe values: {finite}")

    # auto-detect the polarization: TE has Ey only, TM has Ex/Ez instead
    ey = [r[6] for r in fp]
    if max(abs(v) for v in ey) > 0.0:
        sig, pol = ey, "TE (Ey)"
    else:
        sig = [(r[5] ** 2 + r[7] ** 2) ** 0.5 for r in fp]
        pol = "TM (|Ex,Ez|)"
    half = len(sig) // 2
    early = max(abs(v) for v in sig[:half]) if half else 0.0
    late = max(abs(v) for v in sig[half:]) if half else 0.0
    ratio = late / early if early > 0 else float("inf")
    good = ratio <= growth_max
    ok = ok and good
    print(f"{pol} peak: early={early:.4e}  late={late:.4e}  late/early={ratio:.4f} "
          f"({'OK' if good else 'FAIL'}, limit {growth_max})")

    e0 = fe[0][2]
    emax = max(r[2] for r in fe)
    eend = fe[-1][2]
    g = emax / e0 if e0 > 0 else float("inf")
    drift = abs(eend - e0) / e0 if e0 > 0 else float("inf")
    good = (g <= 1.0 + drift_max) and (drift <= drift_max)
    ok = ok and good
    print(f"total energy: E0={e0:.6e}  Emax/E0={g:.6f}  |Eend-E0|/E0={drift:.3e} "
          f"({'OK' if good else 'FAIL'}, limit {drift_max})")

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
