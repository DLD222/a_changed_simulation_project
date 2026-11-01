#!/usr/bin/env python3
"""Kottke 2008 tau-transform unit test (pure Python, no WarpX).

Reference: C. Kottke, A. Farjadpour, S. G. Johnson,
  Phys. Rev. E 77, 036611 (2008), Eq. (4), (22), (23).

Checks:
  1. isotropic reduction -> eps^-1 = P<eps^-1> + (I-P)<eps>^-1
  2. single-material identity -> tau^-1(tau(eps)) = eps (machine precision)
  3. rotation consistency -> rotating materials+normal rotates the result
  4. SPD -> result positive definite for SPD inputs
  5. tau-average property -> tau(eps_bar) = <tau(eps)> (machine precision)
"""
import numpy as np

np.random.seed(20261101)

def frame(n):
    """Orthonormal columns [n, t1, t2]."""
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = a - n * (a @ n); t1 /= np.linalg.norm(t1)
    t2 = np.cross(n, t1)
    return np.stack([n, t1, t2], axis=1)

def tau(eps):
    """Eq. (4). eps in the frame where axis 0 is the interface normal."""
    e = np.asarray(eps, float)
    e11 = e[0, 0]
    T = np.empty((3, 3))
    T[0, 0] = -1.0 / e11
    T[0, 1] = e[0, 1] / e11
    T[0, 2] = e[0, 2] / e11
    T[1, 0] = e[1, 0] / e11
    T[1, 1] = e[1, 1] - e[1, 0] * e[0, 1] / e11
    T[1, 2] = e[1, 2] - e[1, 0] * e[0, 2] / e11
    T[2, 0] = e[2, 0] / e11
    T[2, 1] = e[2, 1] - e[2, 0] * e[0, 1] / e11
    T[2, 2] = e[2, 2] - e[2, 0] * e[0, 2] / e11
    return T

def tau_inv(T):
    """Eq. (23)."""
    T = np.asarray(T, float)
    t11 = T[0, 0]
    e = np.empty((3, 3))
    e[0, 0] = -1.0 / t11
    e[0, 1] = -T[0, 1] / t11
    e[0, 2] = -T[0, 2] / t11
    e[1, 0] = -T[1, 0] / t11
    e[1, 1] = T[1, 1] - T[1, 0] * T[0, 1] / t11
    e[1, 2] = T[1, 2] - T[1, 0] * T[0, 2] / t11
    e[2, 0] = -T[2, 0] / t11
    e[2, 1] = T[2, 1] - T[2, 0] * T[0, 1] / t11
    e[2, 2] = T[2, 2] - T[2, 0] * T[0, 2] / t11
    return e

def subpixel(eps1, eps2, f, n):
    """Effective eps (lab frame) for a two-material cell with fill fraction f
    of material 1 and interface normal n."""
    R = frame(n)
    e1 = R.T @ np.asarray(eps1, float) @ R
    e2 = R.T @ np.asarray(eps2, float) @ R
    Tbar = f * tau(e1) + (1.0 - f) * tau(e2)
    ebar_n = tau_inv(Tbar)
    return R @ ebar_n @ R.T, Tbar, ebar_n

def rand_spd():
    M = np.random.randn(3, 3)
    A = M @ M.T + np.eye(3)
    return A

def rand_rot():
    Q, _ = np.linalg.qr(np.random.randn(3, 3))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q

# ---------------------------------------------------------------- checks
print("=" * 72)
print("Kottke tau-transform unit test")

# 1. isotropic reduction
a, b, f = 1.0, 4.0, 0.3
n = np.array([0.3, -0.5, 0.8]); n = n / np.linalg.norm(n)
eps1 = a * np.eye(3); eps2 = b * np.eye(3)
ebar, Tbar, ebar_n = subpixel(eps1, eps2, f, n)
meps = f * a + (1 - f) * b
minveps = f / a + (1 - f) / b
P = np.outer(n, n)
expected_inv = P * minveps + (np.eye(3) - P) * (1.0 / meps)
err1 = np.linalg.norm(np.linalg.inv(ebar) - expected_inv) / np.linalg.norm(expected_inv)
print(f"[1] isotropic reduction: rel err = {err1:.3e}   (expect ~1e-16)")

# 2. single-material identity
eps = rand_spd(); n = np.array([0.2, 0.9, -0.3]); n /= np.linalg.norm(n)
ebar, _, _ = subpixel(eps, eps.copy(), 0.37, n)
err2 = np.linalg.norm(ebar - eps) / np.linalg.norm(eps)
print(f"[2] single-material identity: rel err = {err2:.3e}   (expect ~1e-15)")

# 3. rotation consistency
eps1 = rand_spd(); eps2 = rand_spd(); f = 0.42
n = np.array([0.1, 0.2, 0.97]); n /= np.linalg.norm(n)
ebar, _, _ = subpixel(eps1, eps2, f, n)
Q = rand_rot()
ebar_rot, _, _ = subpixel(Q @ eps1 @ Q.T, Q @ eps2 @ Q.T, f, Q @ n)
err3 = np.linalg.norm(ebar_rot - Q @ ebar @ Q.T) / np.linalg.norm(ebar)
print(f"[3] rotation consistency: rel err = {err3:.3e}   (expect ~1e-15)")

# 4. SPD
ok = True
for _ in range(200):
    e1, e2 = rand_spd(), rand_spd()
    nn = np.random.randn(3); nn /= np.linalg.norm(nn)
    ff = np.random.rand()
    eb, _, _ = subpixel(e1, e2, ff, nn)
    w = np.linalg.eigvalsh(0.5 * (eb + eb.T))
    if w.min() <= 0 or np.abs(np.linalg.norm(eb - eb.T)) > 1e-12:
        ok = False; break
print(f"[4] SPD (200 random cases): {'PASS' if ok else 'FAIL'}")

# 5. tau-average property: tau(ebar_n) = Tbar
ebar, Tbar, ebar_n = subpixel(rand_spd(), rand_spd(), 0.61,
                              np.array([0.0, 0.6, 0.8]))
err5 = np.linalg.norm(tau(ebar_n) - Tbar) / np.linalg.norm(Tbar)
print(f"[5] tau-average property: rel err = {err5:.3e}   (expect ~1e-15)")

print("=" * 72)
print("PASS" if max(err1, err2, err3, err5) < 1e-12 and ok else "CHECK FAILURES")
