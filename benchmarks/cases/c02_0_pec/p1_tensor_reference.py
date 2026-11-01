#!/usr/bin/env python3
"""P1 tensor reference: independent numpy implementation of

    chi_inv_isotropic / chi_inv_anisotropic (Kottke tau) / chi_bulk / chi_diag

Generates a fixed-seed case set covering L0.5-L0.14 of
Docs/research_notes/P1_BUG_HUNT_PROTOCOL.md, writes p1_cases.csv and
p1_expected.csv, and can compare against p1_cpp.csv.

Usage:
  python3 p1_tensor_reference.py --write          # generate cases + expected
  python3 p1_tensor_reference.py --compare p1_cpp.csv
"""
import argparse
import csv
import numpy as np

np.random.seed(20261101)


# ---------------------------------------------------------------- tensor math
def frame(n):
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = a - n * (a @ n)
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(n, t1)
    return np.stack([n, t1, t2], axis=1)


def tau(e):
    e11 = e[0, 0]
    return np.array([
        [-1.0 / e11, e[0, 1] / e11, e[0, 2] / e11],
        [e[1, 0] / e11, e[1, 1] - e[1, 0] * e[0, 1] / e11,
         e[1, 2] - e[1, 0] * e[0, 2] / e11],
        [e[2, 0] / e11, e[2, 1] - e[2, 0] * e[0, 1] / e11,
         e[2, 2] - e[2, 0] * e[0, 2] / e11],
    ])


def tau_inv(T):
    t11 = T[0, 0]
    return np.array([
        [-1.0 / t11, -T[0, 1] / t11, -T[0, 2] / t11],
        [-T[1, 0] / t11, T[1, 1] - T[1, 0] * T[0, 1] / t11,
         T[1, 2] - T[1, 0] * T[0, 2] / t11],
        [-T[2, 0] / t11, T[2, 1] - T[2, 0] * T[0, 1] / t11,
         T[2, 2] - T[2, 0] * T[0, 2] / t11],
    ])


def _sym(a):
    return 0.5 * (a + a.T)


def chi_inv_anisotropic(e1, e2, f, n):
    """Effective chi^-1 = (tau^-1(f tau(e1')+(1-f) tau(e2')))^-1, in lab frame."""
    R = frame(n)
    e1n = _sym(R.T @ np.asarray(e1, float) @ R)
    e2n = _sym(R.T @ np.asarray(e2, float) @ R)
    T = f * tau(e1n) + (1.0 - f) * tau(e2n)
    ebar_n = tau_inv(T)
    ebar = R @ ebar_n @ R.T
    chi = np.linalg.inv(_sym(ebar))
    return _sym(chi)


def chi_inv_isotropic(eps1, eps2, f, n):
    meps = f * eps1 + (1.0 - f) * eps2
    minveps = f / eps1 + (1.0 - f) / eps2
    P = np.outer(n, n)
    return P * minveps + (np.eye(3) - P) * (1.0 / meps)


def chi_bulk(e):
    return np.linalg.inv(_sym(np.asarray(e, float)))


def chi_diag(eps):
    return np.diag([1.0 / eps] * 3)


# ---------------------------------------------------------------- case set
def rand_spd():
    M = np.random.randn(3, 3)
    return _sym(M @ M.T + np.eye(3))


def rand_rot():
    Q, _ = np.linalg.qr(np.random.randn(3, 3))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


def rand_n(n2d=False):
    if n2d:
        t = np.random.rand() * 2 * np.pi
        return np.array([np.sin(t), 0.0, np.cos(t)])
    v = np.random.randn(3)
    return v / np.linalg.norm(v)


def cases():
    out = []
    f_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
    # (I,I): isotropic degeneracy
    for f in f_vals:
        out.append(("iso_iso", np.eye(3) * 1.0, np.eye(3) * 4.0, f, rand_n()))
        out.append(("iso_iso_2d", np.eye(3) * 1.0, np.eye(3) * 4.0, f, rand_n(True)))
    # (A,I), (I,A), (A,B)
    for f in (0.25, 0.5, 0.75):
        out.append(("aniso_iso", rand_spd(), np.eye(3) * 4.0, f, rand_n()))
        out.append(("iso_aniso", np.eye(3) * 1.0, rand_spd(), f, rand_n()))
        out.append(("aniso_aniso", rand_spd(), rand_spd(), f, rand_n()))
        out.append(("aniso_aniso_2d", rand_spd(), rand_spd(), f, rand_n(True)))
    # rotated diagonal tensors (both sides)
    for f in (0.3, 0.7):
        for _ in range(4):
            Q1, Q2 = rand_rot(), rand_rot()
            d1 = np.diag(np.random.uniform(1.0, 5.0, 3))
            d2 = np.diag(np.random.uniform(1.0, 5.0, 3))
            out.append(("rot_rot", Q1 @ d1 @ Q1.T, Q2 @ d2 @ Q2.T, f, rand_n()))
    # f endpoints with anisotropic materials
    for f in (0.0, 1.0):
        out.append(("endpoint", rand_spd(), rand_spd(), f, rand_n()))
    # extreme contrast
    out.append(("extreme", np.eye(3) * 1.0, np.eye(3) * 1e6, 0.5, rand_n()))
    # n -> -n pair (same materials)
    e1, e2 = rand_spd(), rand_spd()
    nn = rand_n()
    out.append(("nflip", e1, e2, 0.4, nn))
    out.append(("nflip", e1, e2, 0.4, -nn))
    # tangent rotation pair (same n)
    out.append(("tanrot", e1, e2, 0.4, nn))
    out.append(("tanrot", e1, e2, 0.4, nn))
    return out


def mat_to_row(m):
    return [m[i, j] for i in range(3) for j in range(3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--compare")
    args = ap.parse_args()

    if args.write:
        rows = []
        exp = []
        for k, (kind, e1, e2, f, n) in enumerate(cases()):
            chi = chi_inv_anisotropic(e1, e2, f, n)
            rows.append([k, kind, *mat_to_row(_sym(e1)), *mat_to_row(_sym(e2)),
                         *n, f])
            exp.append([k, *mat_to_row(chi)])
        with open("p1_cases.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["case", "kind", "e1_xx", "e1_xy", "e1_xz", "e1_yx", "e1_yy",
                        "e1_yz", "e1_zx", "e1_zy", "e1_zz",
                        "e2_xx", "e2_xy", "e2_xz", "e2_yx", "e2_yy",
                        "e2_yz", "e2_zx", "e2_zy", "e2_zz",
                        "nx", "ny", "nz", "f"])
            w.writerows(rows)
        with open("p1_expected.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["case", *[f"chi_{i}{j}" for i in range(3) for j in range(3)]])
            w.writerows(exp)
        print(f"wrote p1_cases.csv ({len(rows)} cases) and p1_expected.csv")
        # self-checks
        e = rand_spd()
        print("self-check single material:",
              np.max(np.abs(chi_inv_anisotropic(e, e, 0.37, rand_n()) - np.linalg.inv(e))))
        meps = 0.3 * 1.0 + 0.7 * 4.0
        minveps = 0.3 / 1.0 + 0.7 / 4.0
        expected = np.diag([1 / meps, 1 / meps, minveps])
        print("self-check isotropic:",
              np.max(np.abs(chi_inv_isotropic(1.0, 4.0, 0.3, np.array([0.0, 0.0, 1.0]))
                              - expected)))
        return

    if args.compare:
        exp = {}
        with open("p1_expected.csv") as fh:
            for row in csv.DictReader(fh):
                exp[int(row["case"])] = np.array(
                    [float(row[f"chi_{i}{j}"]) for i in range(3) for j in range(3)])
        # tolerance per kind: `extreme` uses a 1e6 contrast, where the closed-form
        # inverse loses ~1e-12 to rounding (the formula is fine; cond(A) ~ 1e6).
        kind_of = {}
        with open("p1_cases.csv") as fh:
            for row in csv.DictReader(fh):
                kind_of[int(row["case"])] = row["kind"]
        tol = {"extreme": 1e-10}
        worst = 0.0
        bad = []
        with open(args.compare) as fh:
            for row in csv.DictReader(fh):
                k = int(row["case"])
                chi = np.array([float(row[f"chi_{i}{j}"]) for i in range(3) for j in range(3)])
                ref = exp[k]
                err = np.linalg.norm(chi - ref) / max(np.linalg.norm(ref), 1e-300)
                worst = max(worst, err)
                if err > tol.get(kind_of.get(k, ""), 1e-12):
                    bad.append((k, err))
        print(f"max relative error = {worst:.3e}; failures = {len(bad)}")
        for k, err in bad[:20]:
            print(f"  case {k}: {err:.3e}")
        return

    print("use --write or --compare")


if __name__ == "__main__":
    main()
