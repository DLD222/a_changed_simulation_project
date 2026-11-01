/* Copyright (c) 2026
 *
 * 10-moment (pressure-tensor / anisotropic Euler) equation object.
 * Physics ported from Gkeyll gkyl_wv_ten_moment (moments/zero/).
 *
 * Conserved: U = [rho, rho*ux, rho*uy, rho*uz,
 *                 Pxx_c, Pxy_c, Pxz_c, Pyy_c, Pyz_c, Pzz_c],
 * where Pij_c = rho*ui*uj + Pij.
 *
 * This is the pure 10-moment (heat-flux-free) hyperbolic flux. The closure
 * (k0 regularization / gradient closure) is applied as a separate parabolic
 * source term, not in this hyperbolic flux.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_EQUATIONS_TENMOMENT_H_
#define MOMENTSOLVER_EQUATIONS_TENMOMENT_H_

#include "Equations/MomentEquation.H"

#include <cmath>
#include <algorithm>

namespace moments {

struct TenMoment
{
    static constexpr int num_eqs = 10;
    static constexpr int num_waves = 5;

    struct Context
    {
        Real k0 = 1.0; //!< heat-flux regularization coefficient (closure step)
    };

    /** Map a symmetric index pair (i,j) to the packed pressure index. */
    MOMENTS_GPU_HOST_DEVICE
    static int p_index (int i, int j) noexcept
    {
        // order: (0,0)=0, (0,1)=1, (0,2)=2, (1,1)=3, (1,2)=4, (2,2)=5
        if (i > j) { const int t = i; i = j; j = t; }
        return (i == 0) ? j : ((i == 1) ? (2 + j) : 5);
    }

    /** Primitive: W = [rho, ux, uy, uz, Pxx, Pxy, Pxz, Pyy, Pyz, Pzz]. */
    MOMENTS_GPU_HOST_DEVICE
    static void conserved_to_primitive (const Real* U, Real* W, const Context& ctx) noexcept
    {
        (void)ctx;
        const Real rho = U[0];
        const Real ux = U[1] / rho;
        const Real uy = U[2] / rho;
        const Real uz = U[3] / rho;
        W[0] = rho;
        W[1] = ux; W[2] = uy; W[3] = uz;
        W[4] = U[4] - rho*ux*ux;             // Pxx
        W[5] = U[5] - rho*ux*uy;             // Pxy
        W[6] = U[6] - rho*ux*uz;             // Pxz
        W[7] = U[7] - rho*uy*uy;             // Pyy
        W[8] = U[8] - rho*uy*uz;             // Pyz
        W[9] = U[9] - rho*uz*uz;             // Pzz
    }

    /** Maximum absolute wave speed in direction `dir` (longitudinal acoustic). */
    MOMENTS_GPU_HOST_DEVICE
    static Real max_wave_speed (const Real* U, int dir, const Context& ctx) noexcept
    {
        const Real rho = U[0];
        const Real u   = U[1 + dir] / rho;
        const Real Pdd = (dir == 0) ? (U[4] - rho*u*u)
                                : ((dir == 1) ? (U[7] - rho*u*u)
                                              : (U[9] - rho*u*u));
        return std::abs(u) + std::sqrt(3.0 * Pdd / rho);
    }

    /**
     * Physical flux in direction `dir`.
     * F_k(Eij) = rho*ui*uj*uk + ui*Pjk + uj*Pik + uk*Pij, with q=0 closure.
     */
    MOMENTS_GPU_HOST_DEVICE
    static void flux (const Real* U, int dir, Real* F, const Context& ctx) noexcept
    {
        (void)ctx;
        const Real rho = U[0];
        const Real u[3] = { U[1]/rho, U[2]/rho, U[3]/rho };
        const Real ud = u[dir];

        // packed pressure tensor (order: xx, xy, xz, yy, yz, zz)
        const Real P[6] = {
            U[4] - rho*u[0]*u[0],  U[5] - rho*u[0]*u[1],  U[6] - rho*u[0]*u[2],
            U[7] - rho*u[1]*u[1],  U[8] - rho*u[1]*u[2],  U[9] - rho*u[2]*u[2]
        };
        // conserved second moments Eij = U[4..9]
        const Real E[6] = { U[4], U[5], U[6], U[7], U[8], U[9] };

        F[0] = rho * ud;
        F[1] = rho * u[0] * ud + P[p_index(0, dir)];
        F[2] = rho * u[1] * ud + P[p_index(1, dir)];
        F[3] = rho * u[2] * ud + P[p_index(2, dir)];

        // component (i,j) ordering for E: 0=xx, 1=xy, 2=xz, 3=yy, 4=yz, 5=zz
        const int ci[6] = {0, 0, 0, 1, 1, 2};
        const int cj[6] = {0, 1, 2, 1, 2, 2};
        for (int n = 0; n < 6; ++n) {
            const int i = ci[n], j = cj[n];
            F[4 + n] = E[n] * ud + u[i] * P[p_index(j, dir)] + u[j] * P[p_index(i, dir)];
        }
    }

    /** Numerical flux: local Lax-Friedrichs (phase-1 default). */
    MOMENTS_GPU_HOST_DEVICE
    static void numerical_flux (const Real* UL, const Real* UR, int dir, Real* F, const Context& ctx) noexcept
    {
        const Real amax = std::max(max_wave_speed(UL, dir, ctx), max_wave_speed(UR, dir, ctx));
        Real fL[num_eqs], fR[num_eqs];
        flux(UL, dir, fL, ctx);
        flux(UR, dir, fR, ctx);
        for (int v = 0; v < num_eqs; ++v) {
            F[v] = 0.5 * (fL[v] + fR[v]) - 0.5 * amax * (UR[v] - UL[v]);
        }
    }
};

} // namespace moments

#endif // MOMENTSOLVER_EQUATIONS_TENMOMENT_H_
