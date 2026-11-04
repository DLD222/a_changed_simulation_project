/* Copyright (c) 2026
 *
 * Explicit 5-moment fluid <-> Maxwell coupling (M2).
 *
 * A single Strang step is:
 *     source(dt/2) -> transport(dt) -> source(dt/2)
 * where "transport" = hyperbolic FV advance of each species (Euler5) and of
 * the field (Maxwell6), and "source" = the pointwise Lorentz force + Ohmic
 * heating on the fluid and the Ampere current on the field (forward Euler).
 *
 * The explicit source is phase-1 (dt must resolve the plasma frequency); a
 * Crank-Nicolson closed-form implicit source is reserved for a later phase.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_COUPLING_MOMENT_COUPLING_H_
#define MOMENTSOLVER_COUPLING_MOMENT_COUPLING_H_

#include "Equations/Euler5Eqn.H"
#include "Equations/MaxwellEqn.H"
#include "FV/MomentFVDriver.H"
#include "Coupling/FluidFriction.H"

#include <vector>
#include <cmath>

namespace moments {

/** Per-species 5-moment parameters. */
struct Species5
{
    Real q;      //!< charge
    Real m;      //!< mass
    Real gamma;  //!< ratio of specific heats
};

/**
 * \brief Explicit coupling source over dt, symplectic Euler.
 *
 * Fluid U = [rho, rho*ux, rho*uy, rho*uz, E_tot] per cell (density is not
 * changed by the source); field = [Ex,Ey,Ez,Bx,By,Bz] per cell.
 *
 * The (momentum, E) subsystem is a harmonic oscillator; a plain forward Euler
 * is unstable for it. We use the symplectic ordering:
 *     E_new = E - dt * J_old / eps0        (Ampere)
 *     u_new = u + dt * (q/m) * (rho E_new + u x B)   (Lorentz, uses E_new)
 * and update the total energy so the internal energy (rho*eps) is conserved,
 * which keeps the pressure non-negative.
 */
inline void apply_coupling_source (Real* field, std::vector<Real*>& U,
                                   const std::vector<Species5>& sp,
                                   int nx, Real dt, Real eps0)
{
    const int ns = static_cast<int>(sp.size());
    for (int i = 0; i < nx; ++i) {
        Real* f = field + 6 * i;
        const Real Bx = f[3], By = f[4], Bz = f[5];

        // 1. old current
        Real Jx = 0.0, Jy = 0.0, Jz = 0.0;
        for (int s = 0; s < ns; ++s) {
            const Real qm = sp[s].q / sp[s].m;
            const Real* u = U[s] + 5 * i;
            Jx += qm * u[1]; Jy += qm * u[2]; Jz += qm * u[3];
        }
        // 2. Ampere
        const Real Ex = f[0] - dt * Jx / eps0;
        const Real Ey = f[1] - dt * Jy / eps0;
        const Real Ez = f[2] - dt * Jz / eps0;
        f[0] = Ex; f[1] = Ey; f[2] = Ez;

        // 3. Lorentz (uses E_new, B unchanged) + 4. consistent energy
        for (int s = 0; s < ns; ++s) {
            Real* u = U[s] + 5 * i;
            const Real rho = u[0];
            const Real qm = sp[s].q / sp[s].m;
            const Real ux0 = u[1], uy0 = u[2], uz0 = u[3];

            u[1] = ux0 + dt * (qm * rho * Ex + qm * (uy0 * Bz - uz0 * By));
            u[2] = uy0 + dt * (qm * rho * Ey + qm * (uz0 * Bx - ux0 * Bz));
            u[3] = uz0 + dt * (qm * rho * Ez + qm * (ux0 * By - uy0 * Bx));

            // kinetic-energy change (keeps rho*eps constant):
            const Real uxa = 0.5 * (ux0 + u[1]);
            const Real uya = 0.5 * (uy0 + u[2]);
            const Real uza = 0.5 * (uz0 + u[3]);
            u[4] += dt * qm * rho * (uxa * Ex + uya * Ey + uza * Ez);
        }
    }
}

/**
 * \brief One Strang time step.
 *
 * \param field  [Ex,Ey,Ez,Bx,By,Bz], nx cells (in/out)
 * \param U      vector of species arrays, each [rho,rho u,E], nx cells (in/out)
 * \param sp     species parameters (q, m, gamma)
 * \param c      speed of light
 */
inline void strang_step (Real* field, std::vector<Real*>& U,
                         const std::vector<Species5>& sp,
                         int nx, Real dx, Real dt, Real eps0, Real c, BC bc,
                         const FrictionParams& friction = FrictionParams{})
{
    // per-species friction data (q, m, EOS for T in the Braginskii model)
    std::vector<Real> q(sp.size()), m(sp.size());
    std::vector<IdealGasEos> ig(sp.size());
    std::vector<Eos> eos(sp.size());
    for (std::size_t s = 0; s < sp.size(); ++s) {
        q[s] = sp[s].q; m[s] = sp[s].m;
        ig[s].gamma = sp[s].gamma;
        eos[s] = make_ideal_gas(&ig[s]);
    }

    apply_coupling_source(field, U, sp, nx, 0.5 * dt, eps0);
    apply_friction(U.data(), static_cast<int>(sp.size()), nx, 0.5 * dt,
                   friction, q.data(), m.data(), eos.data());

    for (std::size_t s = 0; s < sp.size(); ++s) {
        Euler5::Context ectx;
        ectx.gamma = sp[s].gamma;
        FV1D<Euler5>::advance(U[s], nx, dx, dt, ectx, bc);
    }

    Maxwell6::Context mctx;
    mctx.c = c;
    FV1D<Maxwell6>::advance(field, nx, dx, dt, mctx, bc);

    apply_coupling_source(field, U, sp, nx, 0.5 * dt, eps0);
    apply_friction(U.data(), static_cast<int>(sp.size()), nx, 0.5 * dt,
                   friction, q.data(), m.data(), eos.data());
}

} // namespace moments

#endif // MOMENTSOLVER_COUPLING_MOMENT_COUPLING_H_
