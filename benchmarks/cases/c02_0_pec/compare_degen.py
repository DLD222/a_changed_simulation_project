#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Compare two runs of the M2a-1a degeneracy check (ECT vs Yee on a grid-aligned
PEC plane). On such a geometry there are no cut cells, so the conformal ECT
solver must reproduce the staircase Yee solver to round-off.

Usage:  python3 compare_degen.py <ect_dir> <yee_dir> [tol]

Exit status is 0 when every probed quantity agrees within `tol` (relative to the
peak amplitude of that quantity).
"""
import sys
import os

FIELDS = {"Ey": 6, "Ez": 7, "Bx": 8, "By": 9, "Bz": 10}


def load(path):
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith("[") or not line.strip():
                continue
            data.append([float(v) for v in line.split()])
    return data


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    ect_dir, yee_dir = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0e-10

    ok = True
    for probe in ("fp1", "fp2"):
        name = os.path.join("diags", "reducedfiles", probe + ".txt")
        a = load(os.path.join(ect_dir, name))
        b = load(os.path.join(yee_dir, name))
        if len(a) != len(b):
            print(f"{probe}: length mismatch {len(a)} vs {len(b)}")
            ok = False
            continue
        print(f"--- {probe} ({len(a)} samples) ---")
        for fname, col in FIELDS.items():
            peak = max((abs(r[col]) for r in a), default=0.0)
            diff = max((abs(ra[col] - rb[col]) for ra, rb in zip(a, b)), default=0.0)
            rel = diff / peak if peak > 0.0 else diff
            flag = "OK " if rel <= tol else "FAIL"
            print(f"  {flag} {fname:3s} peak={peak:.6e}  max|d|={diff:.6e}  rel={rel:.3e}")
            if rel > tol:
                ok = False

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
