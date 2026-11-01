/* Copyright (c) 2026
 *
 * Ten-moment pressure-tensor closure: k0 hyperbolic regularization
 * (Hakim 2008, Gkeyll-style integrating-factor relaxation).
 *
 * The pure 10-moment (heat-flux-free) system is not hyperbolic; the k0
 * regularization relaxes the anisotropic (trace-free) part of the pressure
 * tensor toward isotropy at rate nu = v_th * k0, while preserving the scalar
 * pressure (trace). Exact exponential (integrating-factor) form:
 *     diagonal:  Pii <- (Pii - p)/exp(nu) + p
 *     off-diag:  Pij <- Pij / exp(nu)
 * with p = (Pxx+Pyy+Pzz)/3, v_th = sqrt(p/rho), exp(nu) = exp(v_th*k0*dt).
 *
 * The gradient closure (q ~ -grad T) is reserved for a later phase.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_COUPLING_TENMOMENT_CLOSURE_H_
#define MOMENTSOLVER_COUPLING_TENMOMENT_CLOSURE_H_

#include "Equations/TenMomentEqn.H"

#include <cmath>

namespace moments {

/** \brief Apply k0 relaxation toward isotropy over dt. No-op when k0 == 0. */
inline void apply_k0_relaxation (Real* U, int nx, Real dt, Real k0)
{
    if (k0 == 0.0) { return; }
    for (int i = 0; i < nx; ++i) {
        Real* u = U + 10 * i;
        const Real rho = u[0];
        const Real mx = u[1], my = u[2], mz = u[3];
        const Real ux = mx / rho, uy = my / rho, uz = mz / rho;

        // thermal pressure tensor
        const Real Pxx = u[4] - rho * ux * ux;
        const Real Pxy = u[5] - rho * ux * uy;
        const Real Pxz = u[6] - rho * ux * uz;
        const Real Pyy = u[7] - rho * uy * uy;
        const Real Pyz = u[8] - rho * uy * uz;
        const Real Pzz = u[9] - rho * uz * uz;

        const Real p_scalar = (Pxx + Pyy + Pzz) / 3.0;
        const Real v_th = std::sqrt(p_scalar / rho);
        const Real inv_exp = 1.0 / std::exp(v_th * k0 * dt);

        const Real nPxx = (Pxx - p_scalar) * inv_exp + p_scalar;
        const Real nPyy = (Pyy - p_scalar) * inv_exp + p_scalar;
        const Real nPzz = (Pzz - p_scalar) * inv_exp + p_scalar;
        const Real nPxy = Pxy * inv_exp;
        const Real nPxz = Pxz * inv_exp;
        const Real nPyz = Pyz * inv_exp;

        u[4] = nPxx + rho * ux * ux;
        u[5] = nPxy + rho * ux * uy;
        u[6] = nPxz + rho * ux * uz;
        u[7] = nPyy + rho * uy * uy;
        u[8] = nPyz + rho * uy * uz;
        u[9] = nPzz + rho * uz * uz;
    }
}

} // namespace moments

#endif // MOMENTSOLVER_COUPLING_TENMOMENT_CLOSURE_H_
