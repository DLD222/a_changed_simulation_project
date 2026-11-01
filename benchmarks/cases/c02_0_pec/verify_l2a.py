#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""L2a: bitwise equivalence of the ECT and Yee solvers on a regular domain.

Both runs use the same setup with the embedded-boundary surface pushed outside
of the domain (`warpx.eb_implicit_function = "z - 100"`), so no cell is cut and
the per-unit-y ECT geometry is uniform. In that limit the conformal Faraday
update has to reproduce the staggered Yee update exactly, i.e. the two runs must
agree to round-off (not merely to discretization error).

Usage:  python3 verify_l2a.py <ect_dir> <yee_dir> [tol]

Exit status is 0 when every compared quantity agrees within `tol`, relative to
the peak amplitude of that quantity.
"""
import os
import sys

FIELD_COLS = {"Ex": 5, "Ey": 6, "Ez": 7, "Bx": 8, "By": 9, "Bz": 10}
ENERGY_COLS = {"total": 2, "E": 3, "B": 4}


def load(path):
    data = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("[") or s.startswith("#"):
                continue
            data.append([float(v) for v in line.split()])
    return data


def compare(name, a, b, cols, tol):
    """Return (ok, lines) for one diagnostic file.

    A component that is identically zero in the exact solution (e.g. Bz for the
    TE polarization used here) only carries a round-off floor, so a purely
    relative criterion is meaningless for it. Such components are judged against
    the amplitude of the other components (``noise`` below), while the physical
    components keep the relative criterion.
    """
    ok = True
    lines = [f"--- {name} ({len(a)} samples) ---"]
    if len(a) != len(b):
        return False, lines + [f"  FAIL length mismatch {len(a)} vs {len(b)}"]
    field_peak = max((abs(r[c]) for r in a for c in cols.values()), default=0.0)
    noise = field_peak * 1.0e-13
    for label, col in cols.items():
        peak = max((abs(r[col]) for r in a), default=0.0)
        diff = max((abs(ra[col] - rb[col]) for ra, rb in zip(a, b)), default=0.0)
        rel = diff / peak if peak > 0.0 else diff
        good = (rel <= tol) or (diff <= noise)
        ok = ok and good
        lines.append(
            f"  {'OK  ' if good else 'FAIL'} {label:5s} peak={peak:.6e} "
            f"max|d|={diff:.6e} rel={rel:.3e}"
        )
    return ok, lines


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    ect_dir, yee_dir = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0e-14

    ok = True
    for probe in ("fp1", "fp2"):
        name = os.path.join("diags", "reducedfiles", probe + ".txt")
        good, lines = compare(
            probe, load(os.path.join(ect_dir, name)), load(os.path.join(yee_dir, name)),
            FIELD_COLS, tol,
        )
        ok = ok and good
        print("\n".join(lines))

    name = os.path.join("diags", "reducedfiles", "fe.txt")
    good, lines = compare(
        "FieldEnergy", load(os.path.join(ect_dir, name)), load(os.path.join(yee_dir, name)),
        ENERGY_COLS, tol,
    )
    ok = ok and good
    print("\n".join(lines))

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
