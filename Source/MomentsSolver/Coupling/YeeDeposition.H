/* Copyright (c) 2026
 *
 * Charge-conserving deposition of a moment fluid onto a Yee grid.
 *
 * Layout (same as WarpX's 2D XZ Yee grid):
 *     rho  : cell centre
 *     Jx   : x-face  (the cell's "lower" x face, i.e. i-1/2)
 *     Jy   : y-face  (j-1/2)
 *     Jz   : z-face  (k-1/2)
 *     div  : (Jhi - Jlo)/h per axis, the same stencil the FV driver uses for rho
 *
 * The point of this file is the *consistency* of J with rho, not the physics of
 * the current.  If J is built as (q/m) * rho * u, the discrete continuity
 * equation only holds to truncation error, because rho in the solver is updated
 * by a flux divergence, not by an advection of rho*u.  Instead we build J from
 * the very same RK-weighted mass flux that updated rho:
 *
 *     U^{n+1} = U^n + dt (1/6 K1 + 1/6 K2 + 2/3 K3),    K = -div F
 *  => (rho^{n+1} - rho^n)/dt + div F = 0                      (exact, algebraic)
 *  => J := (q/m) F  gives  (rho_q^{n+1} - rho_q^n)/dt + div J = 0
 *
 * so a Gauss check on the fluid charge is exact to machine precision rather
 * than approximate.  The flux F is returned by FVND::advance_with_mass_flux.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_COUPLING_YEE_DEPOSITION_H_
#define MOMENTSOLVER_COUPLING_YEE_DEPOSITION_H_

#include "FV/MomentFVNDDriver.H"
#include "Equations/MomentsConfig.H"

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>

namespace moments {

/** Result of one charge-conserving deposition step. */
template <int DIM>
struct YeeDeposit
{
    using Real = moments::Real;

    std::vector<Real> Jlo[3];   //!< current on the cell's lower face, per axis
    std::vector<Real> Jhi[3];   //!< current on the cell's upper face, per axis
    std::vector<Real> rho_n;    //!< charge density before the step (cell centre)
    std::vector<Real> rho_np1;  //!< charge density after the step
};

/**
 * \brief Advance the fluid one step and deposit the charge-conserving current.
 *
 * The state is advanced exactly as FVND::advance would (same RK3 arithmetic),
 * and the current is (q/m) times the RK-weighted mass flux on each face.
 *
 * \return the maximum discrete continuity residual
 *         |(rho^{n+1}-rho^n)/dt + div J| over all cells.
 */
template <typename Eq, int DIM>
inline Real deposit_conservative (FVND<Eq, DIM> const& fv,
                                  Real* U, Real dt,
                                  typename Eq::Context const& ctx,
                                  const BC bc[3],
                                  Real q, Real m,
                                  YeeDeposit<DIM>& out)
{
    using Real = moments::Real;
    static constexpr int N = Eq::num_eqs;
    const std::size_t nc = fv.ncell_total();

    // charge before
    out.rho_n.assign(nc, Real(0));
    for (std::size_t c = 0; c < nc; ++c) { out.rho_n[c] = (q / m) * U[c * N + 0]; }

    // advance, recording the RK-weighted mass flux
    std::vector<Real> Flo[3], Fhi[3];
    fv.advance_with_mass_flux(U, dt, ctx, bc, Flo, Fhi);

    // charge after + current
    out.rho_np1.assign(nc, Real(0));
    for (std::size_t c = 0; c < nc; ++c) { out.rho_np1[c] = (q / m) * U[c * N + 0]; }

    for (int d = 0; d < 3; ++d) {
        if (d < DIM) {
            out.Jlo[d].assign(nc, Real(0));
            out.Jhi[d].assign(nc, Real(0));
            for (std::size_t c = 0; c < nc; ++c) {
                out.Jlo[d][c] = (q / m) * Flo[d][c];
                out.Jhi[d][c] = (q / m) * Fhi[d][c];
            }
        } else {
            out.Jlo[d].assign(nc, Real(0));
            out.Jhi[d].assign(nc, Real(0));
        }
    }

    // continuity residual
    Real worst = Real(0);
    for (std::size_t c = 0; c < nc; ++c) {
        Real res = (out.rho_np1[c] - out.rho_n[c]) / dt;
        for (int d = 0; d < DIM; ++d) {
            const Real h = fv.hcell(d);
            res += (out.Jhi[d][c] - out.Jlo[d][c]) / h;
        }
        worst = std::max(worst, std::abs(res));
    }
    return worst;
}

} // namespace moments

#endif // MOMENTSOLVER_COUPLING_YEE_DEPOSITION_H_
