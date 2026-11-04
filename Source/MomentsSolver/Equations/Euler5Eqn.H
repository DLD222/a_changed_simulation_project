/* Copyright (c) 2026
 *
 * 5-moment (compressible Euler) equation object for the MomentSolver.
 * Physics ported from Gkeyll gkyl_wv_euler (moments/zero/), adapted to
 * AMReX/WarpX conventions.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_EQUATIONS_EULER5_H_
#define MOMENTSOLVER_EQUATIONS_EULER5_H_

#include "Equations/MomentEquation.H"
#include "Equations/EquationOfState.H"

#include <cmath>
#include <algorithm>

namespace moments {

/**
 * \brief 5-moment system: U = [rho, rho*ux, rho*uy, rho*uz, E].
 *
 * Pressure:  p = (gamma-1) * (E - 0.5*|rho u|^2 / rho)
 * Wave speeds (direction d):  u_d, u_d +/- sqrt(gamma p / rho)
 */
struct Euler5
{
    static constexpr int num_eqs = 5;
    static constexpr int num_waves = 3; // contact + two acoustic (entropy/shear degenerate)

    struct Context
    {
        Real gamma = 5.0 / 3.0; //!< ideal-gas ratio of specific heats (fallback)
        Eos  eos;               //!< pluggable EOS (used when use_eos == true)
        bool use_eos = false;   //!< if true, use `eos` instead of the gamma law
    };

    /** Specific internal energy from the conserved state. */
    MOMENTS_GPU_HOST_DEVICE
    static Real internal_energy (const Real* U) noexcept
    {
        const Real rho = U[0];
        const Real u2  = (U[1]*U[1] + U[2]*U[2] + U[3]*U[3]) / (rho*rho);
        return U[4] / rho - 0.5 * u2;
    }

    /** Thermal pressure from the conserved state (ideal gas or pluggable EOS). */
    MOMENTS_GPU_HOST_DEVICE
    static Real pressure (const Real* U, const Context& ctx) noexcept
    {
        const Real rho = U[0];
        if (ctx.use_eos && ctx.eos.valid()) {
            return ctx.eos.p(rho, internal_energy(U));
        }
        const Real u2  = (U[1]*U[1] + U[2]*U[2] + U[3]*U[3]) / (rho*rho);
        return (ctx.gamma - 1.0) * (U[4] - 0.5 * rho * u2);
    }

    /** Maximum absolute wave speed in direction `dir` (for CFL). */
    MOMENTS_GPU_HOST_DEVICE
    static Real max_wave_speed (const Real* U, int dir, const Context& ctx) noexcept
    {
        const Real rho = U[0];
        const Real u   = U[1 + dir] / rho;
        if (ctx.use_eos && ctx.eos.valid()) {
            const Real cs2 = ctx.eos.cs2(rho, internal_energy(U));
            return std::abs(u) + std::sqrt(std::max(cs2, Real(0)));
        }
        const Real p   = pressure(U, ctx);
        return std::abs(u) + std::sqrt(ctx.gamma * p / rho);
    }

    /** Physical flux in direction `dir`. */
    MOMENTS_GPU_HOST_DEVICE
    static void flux (const Real* U, int dir, Real* F, const Context& ctx) noexcept
    {
        const Real rho = U[0];
        const Real ux  = U[1] / rho;
        const Real uy  = U[2] / rho;
        const Real uz  = U[3] / rho;
        const Real p   = pressure(U, ctx);

        const Real u = (dir == 0) ? ux : ((dir == 1) ? uy : uz);

        F[0] = rho * u;
        F[1] = U[1] * u;
        F[2] = U[2] * u;
        F[3] = U[3] * u;
        F[4] = (U[4] + p) * u;
        if (dir == 0) { F[1] += p; }
        if (dir == 1) { F[2] += p; }
        if (dir == 2) { F[3] += p; }
    }

    /** Conserved -> primitive: W = [rho, ux, uy, uz, p]. */
    MOMENTS_GPU_HOST_DEVICE
    static void conserved_to_primitive (const Real* U, Real* W, const Context& ctx) noexcept
    {
        const Real rho = U[0];
        W[0] = rho;
        W[1] = U[1] / rho;
        W[2] = U[2] / rho;
        W[3] = U[3] / rho;
        W[4] = pressure(U, ctx);
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

#endif // MOMENTSOLVER_EQUATIONS_EULER5_H_
