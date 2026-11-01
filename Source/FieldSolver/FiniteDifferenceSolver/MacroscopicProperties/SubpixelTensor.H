/* Copyright 2026
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Kottke 2008 anisotropic subpixel averaging, host/device helpers.
 *
 * The formulas are Eq. (4) and Eq. (23) of
 *   C. Kottke, A. Farjadpour, S. G. Johnson, Phys. Rev. E 77, 036611 (2008)
 * with the tau transform defined in the interface frame (first axis = normal).
 *
 * The isotropic branch uses the closed form
 *   chi^-1_ij = n_i n_j (<1/eps> - 1/<eps>) + delta_ij/<eps>
 * which is algebraically identical to the tau route; it is kept as a fast path so
 * that the isotropic results are bit-for-bit reproducible (P0 regression).
 */
#ifndef WARPX_SUBPIXEL_TENSOR_H_
#define WARPX_SUBPIXEL_TENSOR_H_

#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <cmath>

namespace warpx::subpixel
{

/** Row-major 3x3 matrix (POD, usable inside device lambdas). */
struct Mat3
{
    amrex::Real m[9];

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real& operator() (int r, int c) { return m[3*r + c]; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int r, int c) const { return m[3*r + c]; }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 identity ()
{
    Mat3 a{};
    a.m[0] = a.m[4] = a.m[8] = amrex::Real(1.0);
    return a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 scaled (const Mat3& A, amrex::Real s)
{
    Mat3 B{};
    for (int i = 0; i < 9; ++i) { B.m[i] = A.m[i] * s; }
    return B;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 add (const Mat3& A, const Mat3& B)
{
    Mat3 C{};
    for (int i = 0; i < 9; ++i) { C.m[i] = A.m[i] + B.m[i]; }
    return C;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 mul (const Mat3& A, const Mat3& B)
{
    Mat3 C{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            amrex::Real s = amrex::Real(0.0);
            for (int k = 0; k < 3; ++k) { s += A(i,k) * B(k,j); }
            C(i,j) = s;
        }
    }
    return C;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 transpose (const Mat3& A)
{
    Mat3 B{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) { B(i,j) = A(j,i); }
    }
    return B;
}

/** Rotation whose first column is the interface normal.
 *  In 2D (nx, 0, nz) the second axis is (0,1,0), which is automatically orthogonal,
 *  and the third is n x e2 = (-nz, 0, nx). det = |n|^2 = 1 for a unit normal. */
Mat3 frame (const amrex::Real n_in[3])
{
    // General Gram-Schmidt frame (columns: n, t1, t2), matching the Python
    // reference. The 2D case is just n_y = 0; hard-coding e2 = (0,1,0) broke every
    // 3D case (n_y != 0), which is what p1_tensor_test.cpp caught.
    amrex::Real n[3] = {n_in[0], n_in[1], n_in[2]};
    amrex::Real const nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nn > amrex::Real(0.0)) { n[0] /= nn; n[1] /= nn; n[2] /= nn; }

    amrex::Real a[3] = { (std::abs(n[0]) < amrex::Real(0.9)) ? amrex::Real(1.0) : amrex::Real(0.0),
                         (std::abs(n[0]) < amrex::Real(0.9)) ? amrex::Real(0.0) : amrex::Real(1.0),
                         amrex::Real(0.0) };
    amrex::Real const an = a[0]*n[0] + a[1]*n[1] + a[2]*n[2];
    amrex::Real t1[3] = { a[0] - n[0]*an, a[1] - n[1]*an, a[2] - n[2]*an };
    amrex::Real const t1n = std::sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
    if (t1n > amrex::Real(0.0)) { t1[0] /= t1n; t1[1] /= t1n; t1[2] /= t1n; }
    amrex::Real const t2[3] = { n[1]*t1[2] - n[2]*t1[1],
                                n[2]*t1[0] - n[0]*t1[2],
                                n[0]*t1[1] - n[1]*t1[0] };

    Mat3 R{};
    R(0,0) = n[0];  R(1,0) = n[1];  R(2,0) = n[2];
    R(0,1) = t1[0]; R(1,1) = t1[1]; R(2,1) = t1[2];
    R(0,2) = t2[0]; R(1,2) = t2[1]; R(2,2) = t2[2];
    return R;
}

/** Eq. (4): tau(eps) in the frame where axis 1 is the normal. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 tau_of (const Mat3& e)
{
    amrex::Real const e11 = e(0,0);
    amrex::Real const inv = amrex::Real(1.0) / e11;
    Mat3 T{};
    T(0,0) = -inv;
    T(0,1) =  e(0,1) * inv;
    T(0,2) =  e(0,2) * inv;
    T(1,0) =  e(1,0) * inv;
    T(1,1) =  e(1,1) - e(1,0) * e(0,1) * inv;
    T(1,2) =  e(1,2) - e(1,0) * e(0,2) * inv;
    T(2,0) =  e(2,0) * inv;
    T(2,1) =  e(2,1) - e(2,0) * e(0,1) * inv;
    T(2,2) =  e(2,2) - e(2,0) * e(0,2) * inv;
    return T;
}

/** Eq. (23): the inverse of the tau map (not a matrix inverse).
 *  tau^-1 writes the same structure as tau with the sign of the (1,1) entry flipped:
 *      E_00 = -1/t11            E_0j = -t_0j/t11
 *      E_i0 = -t_i0/t11         E_ij = t_ij - t_i0 t_0j/t11   (i,j >= 1) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 tau_inv_of (const Mat3& T)
{
    amrex::Real const t11 = T(0,0);
    Mat3 E{};
    E(0,0) = amrex::Real(-1.0) / t11;
    E(0,1) = -T(0,1) / t11;
    E(0,2) = -T(0,2) / t11;
    E(1,0) = -T(1,0) / t11;
    E(1,1) = T(1,1) - T(1,0) * T(0,1) / t11;
    E(1,2) = T(1,2) - T(1,0) * T(0,2) / t11;
    E(2,0) = -T(2,0) / t11;
    E(2,1) = T(2,1) - T(2,0) * T(0,1) / t11;
    E(2,2) = T(2,2) - T(2,0) * T(0,2) / t11;
    return E;
}

/** Pointwise isotropy test: diagonal with three equal entries.
 *  Used as the single branch criterion in InitSubpixelRows, so that a scalar material
 *  and an explicitly diagonal tensor take the identical closed-form path (bitwise). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_isotropic (const Mat3& A, amrex::Real& scalar)
{
    bool off_zero = (A(0,1) == amrex::Real(0.0)) && (A(0,2) == amrex::Real(0.0))
                 && (A(1,2) == amrex::Real(0.0));
    bool eq_diag  = (A(0,0) == A(1,1)) && (A(1,1) == A(2,2));
    scalar = A(0,0);
    return off_zero && eq_diag;
}

/** Build a diagonal Mat3 from one scalar. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 diag3 (amrex::Real v)
{
    Mat3 A{};
    A(0,0) = A(1,1) = A(2,2) = v;
    return A;
}

/** Closed-form symmetric 3x3 inverse (cofactor / determinant).
 *  Used by the bulk tensor path; aborts if the determinant is non-positive. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 inv3_sym (const Mat3& A)
{
    amrex::Real const a = A(0,0), b = A(0,1), c = A(0,2);
    amrex::Real const d = A(1,1), e = A(1,2), g = A(2,2);
    amrex::Real const det = a*(d*g - e*e) - b*(b*g - e*c) + c*(b*e - d*c);
    // Sylvester criterion for a symmetric 3x3: a>0, a*d-b*b>0, det>0
    AMREX_ASSERT_WITH_MESSAGE(a > amrex::Real(0.0) && (a*d - b*b) > amrex::Real(0.0)
                              && det > amrex::Real(0.0),
        "subpixel bulk: epsilon/mu tensor is not positive definite (Sylvester)");
    amrex::Real const inv = amrex::Real(1.0) / det;
    Mat3 C{};
    C(0,0) =  (d*g - e*e) * inv;
    C(0,1) = -(b*g - e*c) * inv;
    C(0,2) =  (b*e - d*c) * inv;
    C(1,0) = C(0,1);
    C(1,1) =  (a*g - c*c) * inv;
    C(1,2) = -(a*e - b*c) * inv;
    C(2,0) = C(0,2);
    C(2,1) = C(1,2);
    C(2,2) =  (a*d - b*b) * inv;
    return C;
}

/** Symmetric 3x3 eigen-decomposition by cyclic Jacobi rotations.
 *  `A` must be symmetric; on return `lam` holds the eigenvalues and column `k` of
 *  `V` the corresponding eigenvector (V orthogonal). Used by the sigma exponential,
 *  where every matrix function is f(A) = V diag(f(lam)) V^T. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jacobi_eigh (const Mat3& A_in, amrex::Real lam[3], amrex::Real V[9])
{
    amrex::Real a[3][3];
    amrex::Real sc = amrex::Real(0.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            a[i][j] = amrex::Real(0.5) * (A_in(i,j) + A_in(j,i));
            sc = (std::abs(a[i][j]) > sc) ? std::abs(a[i][j]) : sc;
        }
    }
    // Normalise before iterating: the convergence threshold below is absolute, and
    // eps ~ 1e-11 would otherwise stop at a relative accuracy of ~1e-9.
    if (sc > amrex::Real(0.0)) {
        for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) { a[i][j] /= sc; } }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) { V[3*i + j] = (i == j) ? amrex::Real(1.0) : amrex::Real(0.0); }
    }
    for (int sweep = 0; sweep < 50; ++sweep) {
        amrex::Real const off = a[0][1]*a[0][1] + a[0][2]*a[0][2] + a[1][2]*a[1][2];
        if (off < amrex::Real(1e-31)) { break; }
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(a[p][q]) < amrex::Real(1e-40)) { continue; }
                amrex::Real const theta = (a[q][q] - a[p][p]) / (amrex::Real(2.0) * a[p][q]);
                amrex::Real const t = (theta >= amrex::Real(0.0))
                    ? amrex::Real(1.0) / (theta + std::sqrt(theta*theta + amrex::Real(1.0)))
                    : amrex::Real(-1.0) / (-theta + std::sqrt(theta*theta + amrex::Real(1.0)));
                amrex::Real const c = amrex::Real(1.0) / std::sqrt(t*t + amrex::Real(1.0));
                amrex::Real const sn = t * c;
                amrex::Real const app = a[p][p], aqq = a[q][q], apq = a[p][q];
                a[p][p] = app - t * apq;
                a[q][q] = aqq + t * apq;
                a[p][q] = amrex::Real(0.0);
                a[q][p] = amrex::Real(0.0);
                for (int k = 0; k < 3; ++k) {
                    if (k != p && k != q) {
                        amrex::Real const akp = a[k][p], akq = a[k][q];
                        a[k][p] = c * akp - sn * akq;
                        a[p][k] = a[k][p];
                        a[k][q] = sn * akp + c * akq;
                        a[q][k] = a[k][q];
                    }
                }
                for (int k = 0; k < 3; ++k) {
                    amrex::Real const vkp = V[3*k + p], vkq = V[3*k + q];
                    V[3*k + p] = c * vkp - sn * vkq;
                    V[3*k + q] = sn * vkp + c * vkq;
                }
            }
        }
    }
    lam[0] = a[0][0] * sc; lam[1] = a[1][1] * sc; lam[2] = a[2][2] * sc;
}

/** f(A) for a symmetric A, given its eigen-decomposition. kind: 0 = A^(1/2),
 *  1 = A^(-1/2), 2 = exp(-A dt), 3 = (1-exp(-A dt))/A (-> dt as lambda -> 0). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Mat3 mat3_fun (const amrex::Real lam[3], const amrex::Real V[9], int kind, amrex::Real dt)
{
    amrex::Real g[3];
    for (int i = 0; i < 3; ++i) {
        if (kind == 0) { g[i] = std::sqrt(lam[i]); }
        else if (kind == 1) { g[i] = amrex::Real(1.0) / std::sqrt(lam[i]); }
        else if (kind == 2) { g[i] = std::exp(-lam[i] * dt); }
        else {
            g[i] = (std::abs(lam[i]) > amrex::Real(1e-300))
                 ? (amrex::Real(1.0) - std::exp(-lam[i] * dt)) / lam[i]
                 : dt;
        }
    }
    Mat3 out{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            amrex::Real acc = amrex::Real(0.0);
            for (int k = 0; k < 3; ++k) { acc += V[3*i + k] * g[k] * V[3*j + k]; }
            out(i,j) = acc;
        }
    }
    return out;
}

/** The two static operators of the sigma exponential step.
 *
 *      A  = eps^-1 sigma,  M1 = exp(-A dt),  M2 = (1 - exp(-A dt)) A^-1 eps^-1
 *      E^{n+1} = M1 E^n + M2 curlH
 *
 *  With eps, sigma SPD, B = eps^{-1/2} sigma eps^{-1/2} is symmetric, and both
 *  factors of A become functions of B only. sigma = 0 gives M1 = I, M2 = dt eps^-1,
 *  i.e. the plain Yee update. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void sigma_exponential_matrices (const Mat3& eps, const Mat3& sigma, amrex::Real dt,
                                 Mat3& M1, Mat3& M2)
{
    amrex::Real lam[3], V[9];
    jacobi_eigh(eps, lam, V);
    AMREX_ASSERT_WITH_MESSAGE(lam[0] > amrex::Real(0.0) && lam[1] > amrex::Real(0.0)
                              && lam[2] > amrex::Real(0.0),
        "sigma exponential: epsilon is not positive definite");
    Mat3 const eps_h  = mat3_fun(lam, V, 0, dt);   // eps^{1/2}
    Mat3 const eps_ih = mat3_fun(lam, V, 1, dt);   // eps^{-1/2}

    Mat3 const B = mul(eps_ih, mul(sigma, eps_ih));
    amrex::Real lb[3], Vb[9];
    jacobi_eigh(B, lb, Vb);

    Mat3 const edec = mat3_fun(lb, Vb, 2, dt);     // exp(-B dt)
    Mat3 const eint = mat3_fun(lb, Vb, 3, dt);     // (1-exp(-B dt))/B

    M1 = mul(eps_ih, mul(edec, eps_h));
    M2 = mul(eps_ih, mul(eint, eps_ih));
}

/** Isotropic closed form (identical to the tau route for scalar materials).
 *  chi^-1_ij = n_i n_j (<1/eps> - 1/<eps>) + delta_ij / <eps>  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void chi_inv_isotropic (amrex::Real eps1, amrex::Real eps2, amrex::Real f,
                        const amrex::Real n[3], amrex::Real chi[9])
{
    amrex::Real const e_avg  = f * eps1 + (amrex::Real(1.0) - f) * eps2;
    amrex::Real const ei_avg = f / eps1 + (amrex::Real(1.0) - f) / eps2;
    amrex::Real const diff   = ei_avg - amrex::Real(1.0) / e_avg;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            chi[3*i + j] = n[i] * n[j] * diff
                         + ((i == j) ? amrex::Real(1.0) / e_avg : amrex::Real(0.0));
        }
    }
}

/** Full Kottke subpixel: returns chi^-1 = eps_bar^-1 (lab frame, 9 entries).
 *  `n` must be a unit normal. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void chi_inv_anisotropic (const Mat3& e1, const Mat3& e2, amrex::Real f,
                          const amrex::Real n[3], amrex::Real chi[9])
{
    Mat3 const R  = frame(n);
    Mat3 const Rt = transpose(R);

    // rotate both materials into the interface frame
    Mat3 const e1p = mul(Rt, mul(e1, R));
    Mat3 const e2p = mul(Rt, mul(e2, R));

    // tau-average and invert
    Mat3 const tavg = add(scaled(tau_of(e1p), f),
                          scaled(tau_of(e2p), amrex::Real(1.0) - f));
    Mat3 const ebar_p = tau_inv_of(tavg);

    // Invert in the *interface* frame and rotate back afterwards. R is orthogonal, so
    // this does not change the condition number (cond is rotation invariant); it only
    // saves two rotations and keeps the normal/tangential structure explicit.
    // eps11 = n^T eps n is required by tau; SPD guarantees a positive lower bound but
    // extreme contrasts can make it tiny, so check it explicitly.
    AMREX_ASSERT_WITH_MESSAGE(
        ebar_p(0,0) > amrex::Real(1e-30),
        "subpixel: eps_bar_11 is not positive (the interface-frame normal component "
        "must be positive definite)");

    amrex::Real const a = ebar_p(0,0), b = ebar_p(0,1), c = ebar_p(0,2);
    amrex::Real const d = ebar_p(1,1), e = ebar_p(1,2), g = ebar_p(2,2);
    amrex::Real const det = a*(d*g - e*e) - b*(b*g - e*c) + c*(b*e - d*c);
    AMREX_ASSERT_WITH_MESSAGE(det > amrex::Real(0.0),
        "subpixel: eps_bar is not positive definite in the interface frame");
    amrex::Real const invdet = amrex::Real(1.0) / det;
    Mat3 chi_p{};
    chi_p(0,0) =  (d*g - e*e) * invdet;
    chi_p(0,1) = -(b*g - e*c) * invdet;
    chi_p(0,2) =  (b*e - d*c) * invdet;
    chi_p(1,0) = chi_p(0,1);
    chi_p(1,1) =  (a*g - c*c) * invdet;
    chi_p(1,2) = -(a*e - b*c) * invdet;
    chi_p(2,0) = chi_p(0,2);
    chi_p(2,1) = chi_p(1,2);
    chi_p(2,2) =  (a*d - b*b) * invdet;

    Mat3 const chi_m = mul(R, mul(chi_p, Rt));
    for (int i = 0; i < 9; ++i) { chi[i] = chi_m.m[i]; }
}

} // namespace warpx::subpixel

#endif // WARPX_SUBPIXEL_TENSOR_H_
