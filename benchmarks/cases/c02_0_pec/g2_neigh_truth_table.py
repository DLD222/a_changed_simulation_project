#!/usr/bin/env python3
"""G2 step 1: truth table for the ECT neighbour encoding in 2D XZ.

Motivation
----------
`WarpXFaceExtensions.cpp` stores neighbour offsets as a `Neighbours` enum, which
`FaceInfoBox::uint8_to_inds` decodes back into a pair `(i_n, j_n)`. The pair is
interpreted as offsets *inside the face plane*, so it is dimension dependent:

    dim = 0  (face normal to x):  in-plane directions are (y, z)
    dim = 1  (face normal to y):  in-plane directions are (x, z)
    dim = 2  (face normal to z):  in-plane directions are (x, y)

In the 2D XZ build the physical coordinates are (x, z), stored as (i, j); the k
index is a dummy of extent 1. Therefore only certain combinations are meaningful:

    dim = 0 -> neighbour offset must be (0, j_n)      [only z survives]
    dim = 1 -> neighbour offset is (i_n, j_n)         [x and z both survive]
    dim = 2 -> neighbour offset must be (i_n, 0)      [only x survives]

The current `GetNeigh` has no XZ branch for dim == 0 (it falls back to the 3D
`arr(i, j+i_n, k+j_n)`, which writes into the dummy k direction), and the XZ
branch for dim == 2 uses `arr(i+i_n, j+j_n, k)` (uses j_n, which is a y offset
that does not exist).

This script enumerates all 8 neighbours x 3 normals and prints, for each case,
the offset that is actually accessed today vs the one that is correct, plus a
3x3 grid simulation to show where the neighbour lands.
"""
from itertools import product

# --- encoding used by FaceInfoBox::addConnectedNeighbor -----------------------
ENC = {
    (-1, -1): "nw", (-1, 0): "w", (-1, 1): "sw",
    (0, -1): "n", (0, 1): "s",
    (1, -1): "ne", (1, 0): "e", (1, 1): "se",
}
DEC = {v: k for k, v in ENC.items()}


def current_offset(dim, i_n, j_n):
    """Offset (di, dj) actually accessed by the current GetNeigh in a 2D build."""
    if dim == 0:
        # arr(i, j + i_n, k + j_n): dj = i_n, and j_n goes into the dummy k
        return (0, i_n), ("k+%d" % j_n if j_n else "k",)
    if dim == 1:
        # XZ branch: arr(i + i_n, j + j_n, k)
        return (i_n, j_n), ()
    # dim == 2, XZ branch: arr(i + i_n, j + j_n, k)
    return (i_n, j_n), ()


def expected_offset(dim, i_n, j_n):
    """Offset (di, dj) that keeps the neighbour inside the face plane in XZ."""
    if dim == 0:
        return (0, j_n)
    if dim == 1:
        return (i_n, j_n)
    return (i_n, 0)


def describe(off):
    return f"({off[0]:+d},{off[1]:+d})"


def main():
    print("G2 truth table: ECT neighbours in 2D XZ (i=x, j=z, k=dummy)\n")
    print(f"{'neighbour':10s} {'(i_n,j_n)':10s} {'dim':4s} "
          f"{'current (di,dj)':18s} {'correct (di,dj)':18s} verdict")
    print("-" * 78)

    bad = []
    for (i_n, j_n), name in sorted(ENC.items(), key=lambda kv: kv[0]):
        for dim in (0, 1, 2):
            cur, extra = current_offset(dim, i_n, j_n)
            exp = expected_offset(dim, i_n, j_n)
            ok = (cur == exp) and not extra
            verdict = "OK" if ok else "WRONG"
            if not ok:
                bad.append((name, i_n, j_n, dim, cur, exp, extra))
            print(f"{name:10s} ({i_n:+d},{j_n:+d})  {dim:3d}  "
                  f"{describe(cur):18s} {describe(exp):18s} {verdict}")

    print("-" * 78)
    print(f"{len(bad)} of 24 combinations are wrong\n")

    print("Where the neighbour lands on a 3x3 grid (centre = (1,1)):\n")
    print(f"{'dim':4s} {'neighbour':10s} {'current':12s} {'correct':12s}")
    for dim in (0, 1, 2):
        for (i_n, j_n), name in sorted(ENC.items(), key=lambda kv: kv[0])[:4]:
            cur, _ = current_offset(dim, i_n, j_n)
            exp = expected_offset(dim, i_n, j_n)
            cpos = (1 + cur[0], 1 + cur[1])
            epos = (1 + exp[0], 1 + exp[1])
            print(f"{dim:3d}  {name:10s} {str(cpos):12s} {str(epos):12s}")

    print("\nRound-trip check (encode -> decode must be the identity):")
    for (i, j), name in ENC.items():
        back = DEC[name]
        assert back == (i, j), (i, j, name, back)
    print("  all 8 neighbours round-trip OK")

    print("\nConclusion")
    print("  dim=0: current code writes the j_n offset into the dummy k index.")
    print("         In XZ it must be arr(i, j + j_n, k).")
    print("  dim=2: current code uses j_n (a y offset) as a j offset.")
    print("         In XZ it must be arr(i + i_n, j, k).")
    print("  dim=1: already correct in XZ.")


if __name__ == "__main__":
    main()
