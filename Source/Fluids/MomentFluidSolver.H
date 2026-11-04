/* Copyright (c) 2026
 *
 * Minimal bridge between the MomentSolver 5-moment fluid and WarpX's field
 * solver: deposit a charge-conserving current into `current_fp` so that
 * `EvolveE` consumes it exactly like a particle current.
 *
 * The current is built from the very mass flux that updated the density
 * (see MomentsSolver/Coupling/YeeDeposition.H), so
 *
 *     (rho^{n+1} - rho^n)/dt + div J = 0          (algebraic identity),
 *
 * and because `EvolveE` advances E with that same J at the same staggering,
 *
 *     div E^{n+1} - rho^{n+1}/eps0 = div E^n - rho^n/eps0 .
 *
 * Scope of this prototype (asserted, not assumed):
 *   - 2D XZ, single mesh-refinement level, single box covering the domain;
 *   - periodic in both directions;
 *   - one fluid species, host-side (no GPU);
 *   - the fluid is a *current source only*: it does not feel the fields yet,
 *     so a run measures the deposition path, not a two-way coupling.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MOMENTFLUIDSOLVER_H_
#define WARPX_MOMENTFLUIDSOLVER_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>

/**
 * \brief Adds a 5-moment fluid's charge-conserving current to `current_fp`.
 */
class MomentFluidSolver
{
public:
    MomentFluidSolver () = default;

    /** True if `fluids.moment.enabled = 1` (cached; safe to call every step). */
    static bool IsEnabled ();

    /**
     * \brief Advance the fluid by `dt` and add its current to the Yee faces.
     *
     * \param[in]     dt    time step
     * \param[in,out] J     the three `current_fp` MultiFabs (x-, z-faces)
     * \param[in]     geom  level geometry (single box, periodic)
     * \param[in]     lev   mesh-refinement level (only 0 supported)
     */
    void AdvanceAndDeposit (amrex::Real dt,
                            amrex::MultiFab* J[3],
                            amrex::Geometry const & geom,
                            int lev);

    /** max |div E - rho_q/eps0| on the cell centres (diagnostic). */
    amrex::Real GaussResidual (amrex::MultiFab const & Ex,
                               amrex::MultiFab const & Ez,
                               amrex::Geometry const & geom) const;

    bool Initialized () const { return m_U != nullptr; }

private:
    //! [rho, rho ux, rho uy, rho uz, E] at cell centres, single box
    std::unique_ptr<amrex::MultiFab> m_U;

    // parameters
    amrex::Real m_q     = -1.0;
    amrex::Real m_m     =  1.0;
    amrex::Real m_gamma = 5.0 / 3.0;
    amrex::Real m_eps0  = 1.0;
    amrex::Real m_amp   = 0.01;    //!< density perturbation amplitude
    amrex::Real m_p0    = 1.0e-4;  //!< background pressure

    void InitData (amrex::Geometry const & geom);
};

#endif // WARPX_MOMENTFLUIDSOLVER_H_
