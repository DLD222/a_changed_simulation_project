/* Copyright (c) 2026
 *
 * Vacuum Maxwell equation object (6 variables, linear hyperbolic system).
 *
 * Conserved: U = [Ex, Ey, Ez, Bx, By, Bz].
 * This is the FV equation-object reference for "vacuum EM propagation".
 * The production field solver reuses WarpX's Yee FDTD + divergence cleaning;
 * a divergence-cleaned GLM (8-variable) variant is reserved for later.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_EQUATIONS_MAXWELL_H_
#define MOMENTSOLVER_EQUATIONS_MAXWELL_H_

#include "Equations/MomentEquation.H"

#include <cmath>

namespace moments {

struct Maxwell6
{
    static constexpr int num_eqs = 6;
    static constexpr int num_waves = 2; // +-c (x-direction)

    struct Context
    {
        Real c = 1.0; //!< speed of light in the medium (vacuum: 1/sqrt(eps0*mu0))
    };

    /** Maximum absolute wave speed (isotropic: c). */
    MOMENTS_GPU_HOST_DEVICE
    static Real max_wave_speed (const Real* U, int dir, const Context& ctx) noexcept
    {
        (void)U; (void)dir;
        return ctx.c;
    }

    /**
     * Physical flux in direction `dir` (vacuum Maxwell):
     *   F_x(Ex)=0, F_x(Ey)=+c^2 Bz, F_x(Ez)=-c^2 By,
     *   F_x(Bx)=0, F_x(By)=-Ez,     F_x(Bz)=+Ey.
     */
    MOMENTS_GPU_HOST_DEVICE
    static void flux (const Real* U, int dir, Real* F, const Context& ctx) noexcept
    {
        const Real c2 = ctx.c * ctx.c;
        const Real Ex = U[0], Ey = U[1], Ez = U[2];
        const Real Bx = U[3], By = U[4], Bz = U[5];

        if (dir == 0) { // x-flux:  F(Ey)=+c^2 Bz, F(Ez)=-c^2 By, F(By)=-Ez, F(Bz)=+Ey
            F[0] = 0.0;  F[3] = 0.0;
            F[1] =  c2 * Bz; F[2] = -c2 * By;
            F[4] = -Ez;      F[5] =  Ey;
        } else if (dir == 1) { // y-flux
            F[1] = 0.0;  F[4] = 0.0;
            F[0] = -c2 * Bz; F[2] =  c2 * Bx;
            F[3] =  Ez;      F[5] = -Ex;
        } else { // z-flux
            F[2] = 0.0;  F[5] = 0.0;
            F[0] =  c2 * By; F[1] = -c2 * Bx;
            F[3] = -Ey;      F[4] =  Ex;
        }
    }

    /**
     * Numerical flux: exact upwind (Roe) flux for the linear Maxwell system.
     * |A| is 0 on the longitudinal components (E_dir, B_dir) and c on the
     * transverse components, so no spurious diffusion is added to E_dir/B_dir
     * (important for electrostatic phenomena such as Langmuir oscillations).
     */
    MOMENTS_GPU_HOST_DEVICE
    static void numerical_flux (const Real* UL, const Real* UR, int dir, Real* F, const Context& ctx) noexcept
    {
        const Real c = ctx.c;
        Real fL[num_eqs], fR[num_eqs];
        flux(UL, dir, fL, ctx);
        flux(UR, dir, fR, ctx);
        for (int v = 0; v < num_eqs; ++v) { F[v] = 0.5 * (fL[v] + fR[v]); }
        for (int v = 0; v < num_eqs; ++v) {
            const Real speed = (v == dir || v == 3 + dir) ? Real(0.0) : c;
            F[v] -= 0.5 * speed * (UR[v] - UL[v]);
        }
    }
};

} // namespace moments

#endif // MOMENTSOLVER_EQUATIONS_MAXWELL_H_
