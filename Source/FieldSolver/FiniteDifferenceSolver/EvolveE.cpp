/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <AMReX_Reduce.H>
#include "MacroscopicProperties/MacroscopicProperties.H"
#include <ablastr/coarsen/sample.H>
#include "FiniteDifferenceSolver.H"

#include "Fields.H"
#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER)
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#elif defined(WARPX_DIM_RSPHERE)
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/SphericalYeeAlgorithm.H"
#else
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#endif
#include "EmbeddedBoundary/Enabled.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

using namespace amrex;
using namespace ablastr::fields;

/**
 * \brief Update the E field, over one timestep
 */
void FiniteDifferenceSolver::EvolveE (
    ablastr::fields::MultiFabRegister & fields,
    int lev,
    PatchType patch_type,
    ablastr::fields::VectorField const& Efield,
    std::array< std::unique_ptr<amrex::iMultiFab>,3 > const& eb_update_E,
    amrex::Real const dt
)
{
    using ablastr::fields::Direction;
    using warpx::fields::FieldType;

    const ablastr::fields::VectorField Bfield = patch_type == PatchType::fine ?
        fields.get_alldirs(FieldType::Bfield_fp, lev) : fields.get_alldirs(FieldType::Bfield_cp, lev);
    const ablastr::fields::VectorField Jfield = patch_type == PatchType::fine ?
        fields.get_alldirs(FieldType::current_fp, lev) : fields.get_alldirs(FieldType::current_cp, lev);

    amrex::MultiFab* Ffield = nullptr;
    if (fields.has(FieldType::F_fp, lev)) {
        Ffield = patch_type == PatchType::fine ?
                 fields.get(FieldType::F_fp, lev) : fields.get(FieldType::F_cp, lev);
    }

    // Select algorithm (The choice of algorithm is a runtime option,
    // but we compile code for each algorithm, using templates)
#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER)
    if (m_fdtd_algo == ElectromagneticSolverAlgo::Yee){
        EvolveECylindrical <CylindricalYeeAlgorithm> ( Efield, Bfield, Jfield, eb_update_E, Ffield, lev, dt );
#elif defined(WARPX_DIM_RSPHERE)
    if (m_fdtd_algo == ElectromagneticSolverAlgo::Yee){
        amrex::ignore_unused(eb_update_E);
        EvolveESpherical <SphericalYeeAlgorithm> ( Efield, Bfield, Jfield, Ffield, lev, dt );
#else
    if (m_grid_type == GridType::Collocated) {

        EvolveECartesian <CartesianNodalAlgorithm> ( Efield, Bfield, Jfield, eb_update_E, Ffield, lev, dt );

    } else if (m_fdtd_algo == ElectromagneticSolverAlgo::ECT && WarpX::ect_conformal_E > 0) {

        EvolveECartesianECT ( fields, Efield, Bfield, Jfield, eb_update_E, lev, dt, nullptr );

    } else if (m_fdtd_algo == ElectromagneticSolverAlgo::Yee || m_fdtd_algo == ElectromagneticSolverAlgo::ECT) {

        EvolveECartesian <CartesianYeeAlgorithm> ( Efield, Bfield, Jfield, eb_update_E, Ffield, lev, dt );

    } else if (m_fdtd_algo == ElectromagneticSolverAlgo::CKC) {

        EvolveECartesian <CartesianCKCAlgorithm> ( Efield, Bfield, Jfield, eb_update_E, Ffield, lev, dt );

#endif
    } else {
        WARPX_ABORT_WITH_MESSAGE("EvolveE: Unknown algorithm");
    }

}


#if !defined(WARPX_DIM_RZ) && !defined(WARPX_DIM_RCYLINDER) && !defined(WARPX_DIM_RSPHERE)

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveECartesian (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Jfield,
    std::array< std::unique_ptr<amrex::iMultiFab>,3> const& eb_update_E,
    amrex::MultiFab const* Ffield,
    int lev, amrex::Real const dt ) {

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);
    Real constexpr c2 = PhysConst::c2;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        auto wt = static_cast<amrex::Real>(amrex::second());

        // Extract field data for this grid/tile
        Array4<Real> const& Ex = Efield[0]->array(mfi);
        Array4<Real> const& Ey = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);
        Array4<Real> const& Bx = Bfield[0]->array(mfi);
        Array4<Real> const& By = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);
        Array4<Real> const& jx = Jfield[0]->array(mfi);
        Array4<Real> const& jy = Jfield[1]->array(mfi);
        Array4<Real> const& jz = Jfield[2]->array(mfi);

        // Extract structures indicating whether the E field should be updated
        amrex::Array4<int> update_Ex_arr, update_Ey_arr, update_Ez_arr;
        if (EB::enabled()) {
            update_Ex_arr = eb_update_E[0]->array(mfi);
            update_Ey_arr = eb_update_E[1]->array(mfi);
            update_Ez_arr = eb_update_E[2]->array(mfi);
        }

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        auto const n_coefs_x = static_cast<int>(m_stencil_coefs_x.size());
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        auto const n_coefs_y = static_cast<int>(m_stencil_coefs_y.size());
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        auto const n_coefs_z = static_cast<int>(m_stencil_coefs_z.size());

        // Extract tileboxes for which to loop
        Box const& tex  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tey  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tez  = mfi.tilebox(Efield[2]->ixType().toIntVect());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tex, tey, tez,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                // Skip field push in the embedded boundaries
                if (update_Ex_arr && update_Ex_arr(i, j, k) == 0) { return; }

                Ex(i, j, k) += c2 * dt * (
                    - T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k)
                    + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k)
                    - PhysConst::mu0 * jx(i, j, k) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                // Skip field push in the embedded boundaries
                if (update_Ey_arr && update_Ey_arr(i, j, k) == 0) { return; }

                Ey(i, j, k) += c2 * dt * (
                    - T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k)
                    + T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k)
                    - PhysConst::mu0 * jy(i, j, k) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                // Skip field push in the embedded boundaries
                if (update_Ez_arr && update_Ez_arr(i, j, k) == 0) { return; }

                Ez(i, j, k) += c2 * dt * (
                    - T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k)
                    + T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k)
                    - PhysConst::mu0 * jz(i, j, k) );
            }

        );

        // If F is not a null pointer, further update E using the grad(F) term
        // (hyperbolic correction for errors in charge conservation)
        if (Ffield) {

            // Extract field data for this grid/tile
            const Array4<Real const> F = Ffield->array(mfi);

            // Loop over the cells and update the fields
            amrex::ParallelFor(tex, tey, tez,

                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Ex(i, j, k) += c2 * dt * T_Algo::UpwardDx(F, coefs_x, n_coefs_x, i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Ey(i, j, k) += c2 * dt * T_Algo::UpwardDy(F, coefs_y, n_coefs_y, i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Ez(i, j, k) += c2 * dt * T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, k);
                }

            );

        }

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = static_cast<amrex::Real>(amrex::second()) - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }

}

void FiniteDifferenceSolver::EvolveECartesianECT (
    ablastr::fields::MultiFabRegister& fields,
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Jfield,
    std::array< std::unique_ptr<amrex::iMultiFab>,3 > const& eb_update_E,
    int lev, amrex::Real const dt,
    std::unique_ptr<MacroscopicProperties> const& macroscopic_properties )
{
    using warpx::fields::FieldType;

#if defined(WARPX_DIM_XZ)

    Real constexpr c2 = PhysConst::c2;

    // Denominator floors for the dual lengths (G3): Lx and Lz are also the denominators
    // of Ez and Ex, so a cut line of vanishing length would amplify them. Only the
    // denominators are floored; the Ey flux weights keep the true open lengths.
    auto const cell_size = WarpX::CellSize(lev);
    Real const L_stab_x = static_cast<Real>(WarpX::ect_conformal_E_Astab) * cell_size[0];
    Real const L_stab_z = static_cast<Real>(WarpX::ect_conformal_E_Astab) * cell_size[2];
    auto const dx = cell_size[0];
    auto const dz = cell_size[2];

    // Material support (M2b-1). When macroscopic_properties is null this is exactly
    // vacuum: eps_r = 1 and sigma = 0, i.e. the existing behaviour.
    // vacuum permittivity from mu0 and c2 (PhysConst does not expose ep0 here)
    constexpr amrex::Real epsilon0 = amrex::Real(1.0)/(PhysConst::mu0*PhysConst::c2);
    bool const has_macro = (macroscopic_properties != nullptr);
    // The material MultiFabs are cell-centred; the E components live on staggered
    // positions, so epsilon/sigma must be interpolated exactly as
    // MacroscopicEvolveECartesian does (ablastr::coarsen::sample::Interp).
    amrex::GpuArray<int,3> eps_stag{}, sig_stag{}, mu_stag{}, Ex_stag{}, Ey_stag{}, Ez_stag{}, macro_cr{};
    if (has_macro) {
        eps_stag  = macroscopic_properties->epsilon_IndexType;
        mu_stag   = macroscopic_properties->mu_IndexType;
        sig_stag  = macroscopic_properties->sigma_IndexType;
        Ex_stag   = macroscopic_properties->Ex_IndexType;
        Ey_stag   = macroscopic_properties->Ey_IndexType;
        Ez_stag   = macroscopic_properties->Ez_IndexType;
        macro_cr  = macroscopic_properties->macro_cr_ratio;
    }
    amrex::MultiFab const* eps_mf = has_macro ? &macroscopic_properties->getepsilon_mf() : nullptr;
    amrex::MultiFab const* sigma_mf = has_macro ? &macroscopic_properties->getsigma_mf() : nullptr;
    amrex::MultiFab const* mu_mf = has_macro ? &macroscopic_properties->getmu_mf() : nullptr;

    amrex::MultiFab const& dual_A_y = *fields.get(FieldType::dual_A_y, lev);
    amrex::MultiFab const& dual_L_x = *fields.get(FieldType::dual_L_x, lev);
    amrex::MultiFab const& dual_L_z = *fields.get(FieldType::dual_L_z, lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        Box const& tex = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tey = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tez = mfi.tilebox(Efield[2]->ixType().toIntVect());

        Array4<Real> const& Ex = Efield[0]->array(mfi);
        Array4<Real> const& Ey = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);
        Array4<Real const> const& Bx = Bfield[0]->const_array(mfi);
        Array4<Real const> const& By = Bfield[1]->const_array(mfi);
        Array4<Real const> const& Bz = Bfield[2]->const_array(mfi);
        Array4<Real const> const& jx = Jfield[0]->const_array(mfi);
        Array4<Real const> const& jy = Jfield[1]->const_array(mfi);
        Array4<Real const> const& jz = Jfield[2]->const_array(mfi);

        // Dual geometry, see warpx::embedded_boundary::ComputeDualGeometry
        Array4<Real const> const& Ay = dual_A_y.const_array(mfi);
        Array4<Real const> const& Lx = dual_L_x.const_array(mfi);
        Array4<Real const> const& Lz = dual_L_z.const_array(mfi);

        amrex::Array4<amrex::Real const> eps_arr, sig_arr, mu_arr;
        if (has_macro) {
            eps_arr = eps_mf->const_array(mfi);
            sig_arr = sigma_mf->const_array(mfi);
            mu_arr  = mu_mf->const_array(mfi);
        }

        amrex::Array4<int> update_Ex_arr, update_Ey_arr, update_Ez_arr;
        if (EB::enabled()) {
            update_Ex_arr = eb_update_E[0]->array(mfi);
            update_Ey_arr = eb_update_E[1]->array(mfi);
            update_Ez_arr = eb_update_E[2]->array(mfi);
        }

        // Ex : dual face normal to x, area = Lz (per unit length in y)
        amrex::ParallelFor(tex, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // NOTE: unlike the Yee kernel above (which only returns), the conformal
            // kernel also zeroes masked points, so that "covered" is really E = 0
            // independently of whatever value the field had before.
            if (update_Ex_arr && update_Ex_arr(i, j, k) == 0) { Ex(i, j, k) = amrex::Real(0.0); return; }
            amrex::Real const L = Lz(i, j, k);
            if (L <= amrex::Real(0.0)) { Ex(i, j, k) = amrex::Real(0.0); return; }
            amrex::Real const Ld = (L > L_stab_z) ? L : L_stab_z;
            amrex::Real const dE = -c2 * dt * ( By(i, j, k) - By(i, j-1, k) ) / Ld
                                   - c2 * dt * PhysConst::mu0 * jx(i, j, k);
            if (has_macro && eps_arr(i, j, k) > amrex::Real(0.0)) {
                // M2b-1: E^{n+1} = alpha E^n + beta (curl B / mu0), with the CN
                // coefficients of the macroscopic solver (sigma handled implicitly).
                amrex::Real const eps_abs = amrex::max(ablastr::coarsen::sample::Interp(eps_arr, eps_stag, Ex_stag, macro_cr, i, j, k, 0), amrex::Real(1e-30));
                amrex::Real const sg = amrex::max(amrex::Real(0.0), ablastr::coarsen::sample::Interp(
                    sig_arr, sig_stag, Ex_stag, macro_cr, i, j, k, 0));
                amrex::Real const s_ = sg * dt / (amrex::Real(2.0) * eps_abs);
                amrex::Real const alpha_ = (amrex::Real(1.0) - s_) / (amrex::Real(1.0) + s_);
                amrex::Real const inv_epsr = epsilon0 / eps_abs;
                                amrex::Real const mu_abs = ablastr::coarsen::sample::Interp(
                    mu_arr, mu_stag, Ex_stag, macro_cr, i, j, k, 0);
                amrex::Real const mur_inv = PhysConst::mu0 / amrex::max(mu_abs, amrex::Real(1e-30));
Ex(i, j, k) = alpha_ * Ex(i, j, k) + inv_epsr * mur_inv * dE / (amrex::Real(1.0) + s_);
            } else {
                Ex(i, j, k) += dE;
            }
        });

        // Ez : dual face normal to z, area = Lx (per unit length in y)
        amrex::ParallelFor(tez, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if (update_Ez_arr && update_Ez_arr(i, j, k) == 0) { Ez(i, j, k) = amrex::Real(0.0); return; }
            amrex::Real const L = Lx(i, j, k);
            if (L <= amrex::Real(0.0)) { Ez(i, j, k) = amrex::Real(0.0); return; }
            amrex::Real const Ld = (L > L_stab_x) ? L : L_stab_x;
            amrex::Real const dE = c2 * dt * ( By(i, j, k) - By(i-1, j, k) ) / Ld
                                   - c2 * dt * PhysConst::mu0 * jz(i, j, k);
            if (has_macro && eps_arr(i, j, k) > amrex::Real(0.0)) {
                amrex::Real const eps_abs = amrex::max(ablastr::coarsen::sample::Interp(
                    eps_arr, eps_stag, Ez_stag, macro_cr, i, j, k, 0), amrex::Real(1e-30));
                amrex::Real const sg = amrex::max(amrex::Real(0.0), ablastr::coarsen::sample::Interp(
                    sig_arr, sig_stag, Ez_stag, macro_cr, i, j, k, 0));
                amrex::Real const s_ = sg * dt / (amrex::Real(2.0) * eps_abs);
                amrex::Real const alpha_ = (amrex::Real(1.0) - s_) / (amrex::Real(1.0) + s_);
                amrex::Real const inv_epsr = epsilon0 / eps_abs;
                                amrex::Real const mu_abs = ablastr::coarsen::sample::Interp(
                    mu_arr, mu_stag, Ez_stag, macro_cr, i, j, k, 0);
                amrex::Real const mur_inv = PhysConst::mu0 / amrex::max(mu_abs, amrex::Real(1e-30));
Ez(i, j, k) = alpha_ * Ez(i, j, k) + inv_epsr * mur_inv * dE / (amrex::Real(1.0) + s_);
            } else {
                Ez(i, j, k) += dE;
            }
        });

        // Ey : dual face normal to y, area = Ay (four quadrants of a node)
        amrex::ParallelFor(tey, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if (update_Ey_arr && update_Ey_arr(i, j, k) == 0) { Ey(i, j, k) = amrex::Real(0.0); return; }
            amrex::Real const A = Ay(i, j, k);
            if (A <= amrex::Real(0.0)) { Ey(i, j, k) = amrex::Real(0.0); return; }
            amrex::Real const dE = c2 * dt * (
                  Bx(i, j,   k) * Lx(i, j,   k) - Bx(i, j-1, k) * Lx(i, j-1, k)
                - Bz(i, j,   k) * Lz(i, j,   k) + Bz(i-1, j, k) * Lz(i-1, j, k) ) / A
                - c2 * dt * PhysConst::mu0 * jy(i, j, k);
            if (has_macro && eps_arr(i, j, k) > amrex::Real(0.0)) {
                amrex::Real const eps_abs = amrex::max(ablastr::coarsen::sample::Interp(
                    eps_arr, eps_stag, Ey_stag, macro_cr, i, j, k, 0), amrex::Real(1e-30));
                amrex::Real const sg = amrex::max(amrex::Real(0.0), ablastr::coarsen::sample::Interp(
                    sig_arr, sig_stag, Ey_stag, macro_cr, i, j, k, 0));
                amrex::Real const s_ = sg * dt / (amrex::Real(2.0) * eps_abs);
                amrex::Real const alpha_ = (amrex::Real(1.0) - s_) / (amrex::Real(1.0) + s_);
                amrex::Real const inv_epsr = epsilon0 / eps_abs;
                                amrex::Real const mu_abs = ablastr::coarsen::sample::Interp(
                    mu_arr, mu_stag, Ey_stag, macro_cr, i, j, k, 0);
                amrex::Real const mur_inv = PhysConst::mu0 / amrex::max(mu_abs, amrex::Real(1e-30));
Ey(i, j, k) = alpha_ * Ey(i, j, k) + inv_epsr * mur_inv * dE / (amrex::Real(1.0) + s_);
            } else {
                Ey(i, j, k) += dE;
            }
        });
    }

#else

    amrex::ignore_unused(fields, Efield, Bfield, Jfield, eb_update_E, lev, dt);
    WARPX_ABORT_WITH_MESSAGE("EvolveECartesianECT is only implemented for 2D (XZ) builds");

#endif
}

#elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER)

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveECylindrical (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Jfield,
    std::array< std::unique_ptr<amrex::iMultiFab>,3 > const& eb_update_E,
    amrex::MultiFab const* Ffield,
    int lev, amrex::Real const dt ) {

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        auto wt = static_cast<amrex::Real>(amrex::second());

        // Extract field data for this grid/tile
        Array4<Real> const& Er = Efield[0]->array(mfi);
        Array4<Real> const& Etheta = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);
        Array4<Real> const& Br = Bfield[0]->array(mfi);
        Array4<Real> const& Btheta = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);
        Array4<Real> const& jr = Jfield[0]->array(mfi);
        Array4<Real> const& jtheta = Jfield[1]->array(mfi);
        Array4<Real> const& jz = Jfield[2]->array(mfi);

        // Extract structures indicating whether the E field should be updated
        amrex::Array4<int> update_Er_arr, update_Etheta_arr, update_Ez_arr;
        if (EB::enabled()) {
            update_Er_arr = eb_update_E[0]->array(mfi);
            update_Etheta_arr = eb_update_E[1]->array(mfi);
            update_Ez_arr = eb_update_E[2]->array(mfi);
        }

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_r = m_stencil_coefs_r.dataPtr();
        auto const n_coefs_r = static_cast<int>(m_stencil_coefs_r.size());
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        auto const n_coefs_z = static_cast<int>(m_stencil_coefs_z.size());

        // Extract cylindrical specific parameters
        Real const dr = m_dr;
        int const nmodes = m_nmodes;
        Real const rmin = m_rmin;

        // Extract tileboxes for which to loop
        Box const& ter  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tet  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tez  = mfi.tilebox(Efield[2]->ixType().toIntVect());

        Real constexpr c2 = PhysConst::c2;

        // Loop over the cells and update the fields
        amrex::ParallelFor(ter, tet, tez,

            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){

                // Skip field push in the embedded boundaries
                if (update_Er_arr && update_Er_arr(i, j, 0) == 0) { return; }

                Real const r = rmin + (i + 0.5_rt)*dr; // r on cell-centered point (Er is cell-centered in r)
                Er(i, j, 0, 0) +=  c2 * dt*(
                    - T_Algo::DownwardDz(Btheta, coefs_z, n_coefs_z, i, j, 0, 0)
                    - PhysConst::mu0 * jr(i, j, 0, 0) ); // Mode m=0
                for (int m=1; m<nmodes; m++) { // Higher-order modes
                    Er(i, j, 0, 2*m-1) += c2 * dt*(
                        - T_Algo::DownwardDz(Btheta, coefs_z, n_coefs_z, i, j, 0, 2*m-1)
                        + m * Bz(i, j, 0, 2*m  )/r
                        - PhysConst::mu0 * jr(i, j, 0, 2*m-1) );  // Real part
                    Er(i, j, 0, 2*m  ) += c2 * dt*(
                        - T_Algo::DownwardDz(Btheta, coefs_z, n_coefs_z, i, j, 0, 2*m  )
                        - m * Bz(i, j, 0, 2*m-1)/r
                        - PhysConst::mu0 * jr(i, j, 0, 2*m  ) ); // Imaginary part
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){

                // Skip field push in the embedded boundaries
                if (update_Etheta_arr && update_Etheta_arr(i, j, 0) == 0) { return; }

                Real const r = rmin + i*dr; // r on a nodal grid (Etheta is nodal in r)
                if (r != 0) { // Off-axis, regular Maxwell equations
                    Etheta(i, j, 0, 0) += c2 * dt*(
                        - T_Algo::DownwardDr(Bz, coefs_r, n_coefs_r, i, j, 0, 0)
                        + T_Algo::DownwardDz(Br, coefs_z, n_coefs_z, i, j, 0, 0)
                        - PhysConst::mu0 * jtheta(i, j, 0, 0 ) ); // Mode m=0
                    for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                        Etheta(i, j, 0, 2*m-1) += c2 * dt*(
                            - T_Algo::DownwardDr(Bz, coefs_r, n_coefs_r, i, j, 0, 2*m-1)
                            + T_Algo::DownwardDz(Br, coefs_z, n_coefs_z, i, j, 0, 2*m-1)
                            - PhysConst::mu0 * jtheta(i, j, 0, 2*m-1) ); // Real part
                        Etheta(i, j, 0, 2*m  ) += c2 * dt*(
                            - T_Algo::DownwardDr(Bz, coefs_r, n_coefs_r, i, j, 0, 2*m  )
                            + T_Algo::DownwardDz(Br, coefs_z, n_coefs_z, i, j, 0, 2*m  )
                            - PhysConst::mu0 * jtheta(i, j, 0, 2*m  ) ); // Imaginary part
                    }
                } else { // r==0: on-axis corrections
                    // Ensure that Etheta remains 0 on axis (except for m=1)
                    Etheta(i, j, 0, 0) = 0.; // Mode m=0
                    for (int m=1; m<nmodes; m++) { // Higher-order modes
                        if (m == 1){
                            // The bulk equation could in principle be used here since it does not diverge
                            // on axis. However, it typically gives poor results e.g. for the propagation
                            // of a laser pulse (the field is spuriously reduced on axis). For this reason
                            // a modified on-axis condition is used here: we use the fact that
                            // Etheta(r=0,m=1) should equal -iEr(r=0,m=1), for the fields Er and Etheta to be
                            // independent of theta at r=0. Now with linear interpolation:
                            // Er(r=0,m=1) = 0.5*[Er(r=dr/2,m=1) + Er(r=-dr/2,m=1)]
                            // And using the rule applying for the guards cells
                            // Er(r=-dr/2,m=1) = Er(r=dr/2,m=1). Thus: Etheta(i,j,m) = -i*Er(i,j,m)
                            Etheta(i,j,0,2*m-1) =  Er(i,j,0,2*m  );
                            Etheta(i,j,0,2*m  ) = -Er(i,j,0,2*m-1);
                        } else {
                            Etheta(i, j, 0, 2*m-1) = 0.;
                            Etheta(i, j, 0, 2*m  ) = 0.;
                        }
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){

                // Skip field push in the embedded boundaries
                if (update_Ez_arr && update_Ez_arr(i, j, 0) == 0) { return; }

                Real const r = rmin + i*dr; // r on a nodal grid (Ez is nodal in r)
                if (r != 0) { // Off-axis, regular Maxwell equations
                    Ez(i, j, 0, 0) += c2 * dt*(
                       T_Algo::DownwardDrr_over_r(Btheta, r, dr, coefs_r, n_coefs_r, i, j, 0, 0)
                        - PhysConst::mu0 * jz(i, j, 0, 0  ) ); // Mode m=0
                    for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                        Ez(i, j, 0, 2*m-1) += c2 * dt *(
                            - m * Br(i, j, 0, 2*m  )/r
                            + T_Algo::DownwardDrr_over_r(Btheta, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m-1)
                            - PhysConst::mu0 * jz(i, j, 0, 2*m-1) ); // Real part
                        Ez(i, j, 0, 2*m  ) += c2 * dt *(
                            m * Br(i, j, 0, 2*m-1)/r
                            + T_Algo::DownwardDrr_over_r(Btheta, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m  )
                            - PhysConst::mu0 * jz(i, j, 0, 2*m  ) ); // Imaginary part
                    }
                } else { // r==0: on-axis corrections
                    // For m==0, Btheta is linear in r, for small r
                    // Therefore, the formula below regularizes the singularity
                    Ez(i, j, 0, 0) += c2 * dt*(
                         4*Btheta(i, j, 0, 0)/dr // regularization
                         - PhysConst::mu0 * jz(i, j, 0, 0  ) );
                    // Ensure that Ez remains 0 for higher-order modes
                    for (int m=1; m<nmodes; m++) {
                        Ez(i, j, 0, 2*m-1) = 0.;
                        Ez(i, j, 0, 2*m  ) = 0.;
                    }
                }
            }

        ); // end of loop over cells

        // If F is not a null pointer, further update E using the grad(F) term
        // (hyperbolic correction for errors in charge conservation)
        if (Ffield) {

            // Extract field data for this grid/tile
            const Array4<Real const> F = Ffield->array(mfi);

            // Loop over the cells and update the fields
            amrex::ParallelFor(ter, tet, tez,

                [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){
                    Er(i, j, 0, 0) += c2 * dt * T_Algo::UpwardDr(F, coefs_r, n_coefs_r, i, j, 0, 0);
                    for (int m=1; m<nmodes; m++) { // Higher-order modes
                        Er(i, j, 0, 2*m-1) += c2 * dt * T_Algo::UpwardDr(F, coefs_r, n_coefs_r, i, j, 0, 2*m-1); // Real part
                        Er(i, j, 0, 2*m  ) += c2 * dt * T_Algo::UpwardDr(F, coefs_r, n_coefs_r, i, j, 0, 2*m  ); // Imaginary part
                    }
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){
                    // Mode m=0: no update
                    Real const r = rmin + i*dr; // r on a nodal grid (Etheta is nodal in r)
                    if (r != 0){ // Off-axis, regular Maxwell equations
                        for (int m=1; m<nmodes; m++) { // Higher-order modes
                            Etheta(i, j, 0, 2*m-1) += c2 * dt *  m * F(i, j, 0, 2*m  )/r; // Real part
                            Etheta(i, j, 0, 2*m  ) += c2 * dt * -m * F(i, j, 0, 2*m-1)/r; // Imaginary part
                        }
                    } else { // r==0: on-axis corrections
                        // For m==1, F is linear in r, for small r
                        // Therefore, the formula below regularizes the singularity
                        if (nmodes >= 2) { // needs to have at least m=0 and m=1
                            int const m=1;
                            Etheta(i, j, 0, 2*m-1) += c2 * dt *  m * F(i+1, j, 0, 2*m  )/dr; // Real part
                            Etheta(i, j, 0, 2*m  ) += c2 * dt * -m * F(i+1, j, 0, 2*m-1)/dr; // Imaginary part
                        }
                    }
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){
                    Ez(i, j, 0, 0) += c2 * dt * T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, 0, 0);
                    for (int m=1; m<nmodes; m++) { // Higher-order modes
                        Ez(i, j, 0, 2*m-1) += c2 * dt * T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, 0, 2*m-1); // Real part
                        Ez(i, j, 0, 2*m  ) += c2 * dt * T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, 0, 2*m  ); // Imaginary part
                    }
                }

            ); // end of loop over cells

        } // end of if condition for F

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = static_cast<amrex::Real>(amrex::second()) - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    } // end of loop over grid/tiles

}

#elif defined(WARPX_DIM_RSPHERE)

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveESpherical (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Jfield,
    amrex::MultiFab const* Ffield,
    int lev, amrex::Real const dt ) {

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        auto wt = static_cast<amrex::Real>(amrex::second());

        // Extract field data for this grid/tile
        Array4<Real> const& Er = Efield[0]->array(mfi);
        Array4<Real> const& Etheta = Efield[1]->array(mfi);
        Array4<Real> const& Ephi = Efield[2]->array(mfi);
        Array4<Real> const& Btheta = Bfield[1]->array(mfi);
        Array4<Real> const& Bphi = Bfield[2]->array(mfi);
        Array4<Real> const& jr = Jfield[0]->array(mfi);
        Array4<Real> const& jtheta = Jfield[1]->array(mfi);
        Array4<Real> const& jphi = Jfield[2]->array(mfi);

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_r = m_stencil_coefs_r.dataPtr();
        auto const n_coefs_r = static_cast<int>(m_stencil_coefs_r.size());

        // Extract spherical specific parameters
        Real const dr = m_dr;
        Real const rmin = m_rmin;

        // Extract tileboxes for which to loop
        Box const& ter  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tet  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tep  = mfi.tilebox(Efield[2]->ixType().toIntVect());

        Real constexpr c2 = PhysConst::c2;

        // Loop over the cells and update the fields
        amrex::ParallelFor(ter, tet, tep,

            [=] AMREX_GPU_DEVICE (int i, int /*j*/, int /*k*/){
                Er(i, 0, 0, 0) += c2 * dt*(
                    - PhysConst::mu0 * jr(i, 0, 0, 0) );
            },

            [=] AMREX_GPU_DEVICE (int i, int /*j*/, int /*k*/){
                Real const r = rmin + i*dr; // r on a nodal grid (Etheta is nodal in r)
                if (r != 0) { // Off-axis, regular Maxwell equations
                    Etheta(i, 0, 0, 0) += c2 * dt*(
                        - T_Algo::DownwardDrr_over_r(Bphi, r, dr, coefs_r, n_coefs_r, i, 0, 0, 0)
                        - PhysConst::mu0 * jtheta(i, 0, 0, 0 ) );
                } else { // r==0: on-axis corrections
                    // Ensure that Etheta remains 0 on axis
                    Etheta(i, 0, 0, 0) = 0.;
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int /*j*/, int /*k*/){
                Real const r = rmin + i*dr; // r on a nodal grid (Ephi is nodal in r)
                if (r != 0) { // Off-axis, regular Maxwell equations
                    Ephi(i, 0, 0, 0) += c2 * dt*(
                        + T_Algo::DownwardDrr_over_r(Btheta, r, dr, coefs_r, n_coefs_r, i, 0, 0, 0)
                        - PhysConst::mu0 * jphi(i, 0, 0, 0  ) );
                } else { // r==0: on-axis corrections
                    // Ensure that Ephi remains 0 on axis
                    Ephi(i, 0, 0, 0) = 0.;
                }
            }

        ); // end of loop over cells

        // If F is not a null pointer, further update E using the grad(F) term
        // (hyperbolic correction for errors in charge conservation)
        if (Ffield) {

            // Extract field data for this grid/tile
            const Array4<Real const> F = Ffield->array(mfi);

            // Loop over the cells and update the fields
            amrex::ParallelFor(ter,

                [=] AMREX_GPU_DEVICE (int i, int /*j*/, int /*k*/){
                    Er(i, 0, 0, 0) += c2 * dt * T_Algo::UpwardDr(F, coefs_r, n_coefs_r, i, 0, 0, 0);
                }

            ); // end of loop over cells

        } // end of if condition for F

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = static_cast<amrex::Real>(amrex::second()) - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    } // end of loop over grid/tiles

}


#endif // corresponds to ifndef WARPX_DIM_RZ
