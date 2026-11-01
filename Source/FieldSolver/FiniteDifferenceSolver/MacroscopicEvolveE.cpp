#include <cmath>
#include "FiniteDifferenceSolver.H"

#if !defined(WARPX_DIM_RZ) && !defined(WARPX_DIM_RCYLINDER) && !defined(WARPX_DIM_RSPHERE)
    // currently works only for 3D
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/FieldAccessorFunctors.H"
#endif
#include "EmbeddedBoundary/Enabled.H"
#include "MacroscopicProperties/MacroscopicProperties.H"
#include "WarpX.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"

#include <ablastr/coarsen/sample.H>

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

using namespace amrex;
using namespace ablastr::fields;

void FiniteDifferenceSolver::MacroscopicEvolveE (
    ablastr::fields::MultiFabRegister& fields,
    const MacroscopicSolverAlgo macroscopic_solver_algo,
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Jfield,
    std::array< std::unique_ptr<amrex::iMultiFab>,3 > const& eb_update_E,
    amrex::Real const dt,
    std::unique_ptr<MacroscopicProperties> const& macroscopic_properties)
{

    // Select algorithm (The choice of algorithm is a runtime option,
    // but we compile code for each algorithm, using templates)
#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    amrex::ignore_unused(macroscopic_solver_algo, Efield, Bfield, Jfield, eb_update_E, dt, macroscopic_properties);

    WARPX_ABORT_WITH_MESSAGE("currently macro E-push does not work for RZ");
#else
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_grid_type != GridType::Collocated, "Macroscopic E field solver does not work on collocated grids");


    if (m_fdtd_algo == ElectromagneticSolverAlgo::ECT && WarpX::ect_conformal_E > 0) {
        EvolveECartesianECT ( fields, Efield, Bfield, Jfield, eb_update_E, 0, dt,
                              macroscopic_properties );
    } else if (m_fdtd_algo == ElectromagneticSolverAlgo::Yee) {

        if (macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff) {

            MacroscopicEvolveECartesian <CartesianYeeAlgorithm, LaxWendroffAlgo>
                       ( Efield, Bfield, Jfield, eb_update_E, dt, macroscopic_properties);

        }
        if (macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler) {

            MacroscopicEvolveECartesian <CartesianYeeAlgorithm, BackwardEulerAlgo>
                       ( Efield, Bfield, Jfield, eb_update_E, dt, macroscopic_properties);

        }

    } else if (m_fdtd_algo == ElectromagneticSolverAlgo::CKC) {

        // Note : EvolveE is the same for CKC and Yee.
        // In the templated Yee and CKC calls, the core operations for EvolveE is the same.
        if (macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff) {

            MacroscopicEvolveECartesian <CartesianCKCAlgorithm, LaxWendroffAlgo>
                       ( Efield, Bfield, Jfield, eb_update_E, dt, macroscopic_properties);

        } else if (macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler) {

            MacroscopicEvolveECartesian <CartesianCKCAlgorithm, BackwardEulerAlgo>
                       ( Efield, Bfield, Jfield, eb_update_E, dt, macroscopic_properties);

        }

    } else {
        WARPX_ABORT_WITH_MESSAGE(
            "MacroscopicEvolveE: Unknown algorithm");
    }
#endif

}


#if !defined(WARPX_DIM_RZ) && !defined(WARPX_DIM_RCYLINDER) && !defined(WARPX_DIM_RSPHERE)

template<typename T_Algo, typename T_MacroAlgo>
void FiniteDifferenceSolver::MacroscopicEvolveECartesian (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Jfield,
    std::array< std::unique_ptr<amrex::iMultiFab>,3 > const& eb_update_E,
    amrex::Real const dt,
    std::unique_ptr<MacroscopicProperties> const& macroscopic_properties)
{
    amrex::MultiFab& sigma_mf = macroscopic_properties->getsigma_mf();
    amrex::MultiFab& epsilon_mf = macroscopic_properties->getepsilon_mf();
    amrex::MultiFab& mu_mf = macroscopic_properties->getmu_mf();

    // Index type required for calling ablastr::coarsen::sample::Interp to interpolate macroscopic
    // properties from their respective staggering to the Ex, Ey, Ez locations
    amrex::GpuArray<int, 3> const& sigma_stag = macroscopic_properties->sigma_IndexType;
    amrex::GpuArray<int, 3> const& epsilon_stag = macroscopic_properties->epsilon_IndexType;
    amrex::GpuArray<int, 3> const& macro_cr     = macroscopic_properties->macro_cr_ratio;
    amrex::GpuArray<int, 3> const& Ex_stag = macroscopic_properties->Ex_IndexType;
    amrex::GpuArray<int, 3> const& Ey_stag = macroscopic_properties->Ey_IndexType;
    amrex::GpuArray<int, 3> const& Ez_stag = macroscopic_properties->Ez_IndexType;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

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

        amrex::Array4<int> update_Ex_arr, update_Ey_arr, update_Ez_arr;
        if (EB::enabled()) {
            update_Ex_arr = eb_update_E[0]->array(mfi);
            update_Ey_arr = eb_update_E[1]->array(mfi);
            update_Ez_arr = eb_update_E[2]->array(mfi);
        }

        // material prop //
        amrex::Array4<amrex::Real> const& sigma_arr = sigma_mf.array(mfi);
        amrex::Array4<amrex::Real> const& eps_arr = epsilon_mf.array(mfi);
        amrex::Array4<amrex::Real> const& mu_arr = mu_mf.array(mfi);

        // P0: per-component eps^-1 rows (3 comps each). When the subpixel path is off
        // we fall back to the scalar 1/eps on the diagonal, i.e. the original behaviour.
        bool const has_sub = macroscopic_properties->has_subpixel_rows();
        amrex::Array4<amrex::Real const> eix, eiy, eiz;
        if (has_sub) {
            eix = macroscopic_properties->geteps_inv_x_mf().const_array(mfi);
            eiy = macroscopic_properties->geteps_inv_y_mf().const_array(mfi);
            eiz = macroscopic_properties->geteps_inv_z_mf().const_array(mfi);
        }

        // P2: sigma > 0 cells take the exact exponential step instead of the
        // backward-Euler alpha/beta update. M1/M2 are cell-centred (collocated),
        // so the three E components read the same matrix at the same index.
        bool const has_sig_mat = macroscopic_properties->has_sigma_matrices();
        if (has_sig_mat) {
            // M1/M2 were built with the dt recorded at initialisation; a different dt
            // here would silently use stale matrices (see 技术路径 §7.7).
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                macroscopic_properties->get_sigma_dt() == dt,
                "macroscopic sigma exponential: dt changed after M1/M2 were built");
        }
        amrex::Array4<amrex::Real const> m1a, m2a;
        if (has_sig_mat) {
            m1a = macroscopic_properties->get_sigma_m1_mf().const_array(mfi);
            m2a = macroscopic_properties->get_sigma_m2_mf().const_array(mfi);
        }
        // Node-centred material data (§4.1): the node holds the ONE value shared by both
        // grid copies. BE must use it for BOTH the skip test and alpha/beta -- if the
        // damping still came from the cell interpolation, an interface node would see
        // ~sigma/2 and the over-damped lock would return in another form.
        amrex::Array4<amrex::Real const> sn0{}, sn1{}, sn2{};
        if (has_sig_mat) {
            sn0 = macroscopic_properties->get_sigma_node_mf(0).const_array(mfi);
            sn1 = macroscopic_properties->get_sigma_node_mf(1).const_array(mfi);
            sn2 = macroscopic_properties->get_sigma_node_mf(2).const_array(mfi);
        }

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        auto const n_coefs_x = static_cast<int>(m_stencil_coefs_x.size());
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        auto const n_coefs_y = static_cast<int>(m_stencil_coefs_y.size());
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        auto const n_coefs_z = static_cast<int>(m_stencil_coefs_z.size());

        // This functor computes Hx = Bx/mu
        // Note that mu is cell-centered here and will be interpolated/averaged
        // to the location where the B-field and H-field are defined
        FieldAccessorMacroscopic const Hx(Bx, mu_arr);
        FieldAccessorMacroscopic const Hy(By, mu_arr);
        FieldAccessorMacroscopic const Hz(Bz, mu_arr);

        // Extract tileboxes for which to loop
        Box const& tex  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tey  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tez  = mfi.tilebox(Efield[2]->ixType().toIntVect());
        // starting component to interpolate macro properties to Ex, Ey, Ez locations
        const int scomp = 0;
        // Loop over the cells and update the fields
        // Ex update
        amrex::ParallelFor(tex, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // Skip field push in the embedded boundaries
            if (update_Ex_arr && update_Ex_arr(i, j, k) == 0) { return; }
            // sigma > 0 cells are handled by the exponential step below; leaving them
            // untouched here keeps E^n available for the matrix product.
            // node sigma drives both the skip and the damping (§4.1)
            amrex::Real const sigma_node = has_sig_mat ? sn0(i, j, k, 0) : amrex::Real(0.0);
            if (has_sig_mat && sigma_node > amrex::Real(0.0)) { return; }   // exponential step handles it

            // Interpolate conductivity, sigma, to Ex position on the grid
            // sigma-free cells must not inherit damping from a neighbouring conductor
            // through the face interpolation: that made the vacuum side of a sigma
            // interface decay and locked E (and B) to zero inside the conductor
            // (c03-E6, 技术路径 §7.11). With sigma = 0 everywhere this is a no-op.
            amrex::Real const sigma_interp = sigma_node;   // node sigma, not a cell interpolation
            // Interpolated permittivity, epsilon, to Ex position on the grid
            amrex::Real const epsilon_interp = ablastr::coarsen::sample::Interp(eps_arr, epsilon_stag,
                                                                                Ex_stag, macro_cr, i, j, k, scomp);
            const amrex::Real alpha = T_MacroAlgo::alpha( sigma_interp, epsilon_interp, dt);
            const amrex::Real beta = T_MacroAlgo::beta( sigma_interp, epsilon_interp, dt);
            // P1: full matrix-vector product E_c = alpha E_c + beta0 sum_d (eps^-1)_{c,d} D_d,
            // beta0 = dt/(1+s) = beta*eps. The three curl components are evaluated at the
            // same index (Meep collocated approximation, O(h) at an interface); J stays on
            // the diagonal because its off-diagonal coupling is second order.
            amrex::Real const Dx = - T_Algo::DownwardDz(Hy, coefs_z, n_coefs_z, i, j, k,0)
                                    + T_Algo::DownwardDy(Hz, coefs_y, n_coefs_y, i, j, k,0);
            amrex::Real const Dy =   T_Algo::DownwardDz(Hx, coefs_z, n_coefs_z, i, j, k,0)
                                    - T_Algo::DownwardDx(Hz, coefs_x, n_coefs_x, i, j, k,0);
            amrex::Real const Dz =   T_Algo::DownwardDx(Hy, coefs_x, n_coefs_x, i, j, k,0)
                                    - T_Algo::DownwardDy(Hx, coefs_y, n_coefs_y, i, j, k,0);
            amrex::Real const beta0 = beta * epsilon_interp;
            amrex::Real const ix0 = has_sub ? eix(i,j,k,0) : amrex::Real(1.0)/epsilon_interp;
            amrex::Real const ix1 = has_sub ? eix(i,j,k,1) : amrex::Real(0.0);
            amrex::Real const ix2 = has_sub ? eix(i,j,k,2) : amrex::Real(0.0);
            Ex(i, j, k) = alpha * Ex(i, j, k)
                        + beta0 * ( ix0*Dx + ix1*Dy + ix2*Dz )
                        - beta0 * ix0 * jx(i, j, k);
        });

        // Ey update
        amrex::ParallelFor(tey, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // Skip field push in the embedded boundaries
            if (update_Ey_arr && update_Ey_arr(i, j, k) == 0) { return; }
            // sigma > 0 cells are handled by the exponential step below; leaving them
            // untouched here keeps E^n available for the matrix product.
            // node sigma drives both the skip and the damping (§4.1)
            amrex::Real const sigma_node = has_sig_mat ? sn1(i, j, k, 0) : amrex::Real(0.0);
            if (has_sig_mat && sigma_node > amrex::Real(0.0)) { return; }   // exponential step handles it

            // Interpolate conductivity, sigma, to Ey position on the grid
            // sigma-free cells must not inherit damping from a neighbouring conductor
            // through the face interpolation: that made the vacuum side of a sigma
            // interface decay and locked E (and B) to zero inside the conductor
            // (c03-E6, 技术路径 §7.11). With sigma = 0 everywhere this is a no-op.
            amrex::Real const sigma_interp = sigma_node;   // node sigma, not a cell interpolation
            // Interpolated permittivity, epsilon, to Ey position on the grid
            amrex::Real const epsilon_interp = ablastr::coarsen::sample::Interp(eps_arr, epsilon_stag,
                                                                                Ey_stag, macro_cr, i, j, k, scomp);
            const amrex::Real alpha = T_MacroAlgo::alpha( sigma_interp, epsilon_interp, dt);
            const amrex::Real beta = T_MacroAlgo::beta( sigma_interp, epsilon_interp, dt);

            amrex::Real const Dy0 =   T_Algo::DownwardDz(Hx, coefs_z, n_coefs_z, i, j, k,0)
                                    - T_Algo::DownwardDx(Hz, coefs_x, n_coefs_x, i, j, k,0);
            amrex::Real const Dx0 = - T_Algo::DownwardDz(Hy, coefs_z, n_coefs_z, i, j, k,0)
                                    + T_Algo::DownwardDy(Hz, coefs_y, n_coefs_y, i, j, k,0);
            amrex::Real const Dz0 =   T_Algo::DownwardDx(Hy, coefs_x, n_coefs_x, i, j, k,0)
                                    - T_Algo::DownwardDy(Hx, coefs_y, n_coefs_y, i, j, k,0);
            amrex::Real const beta0 = beta * epsilon_interp;
            amrex::Real const iy0 = has_sub ? eiy(i,j,k,0) : amrex::Real(0.0);
            amrex::Real const iy1 = has_sub ? eiy(i,j,k,1) : amrex::Real(1.0)/epsilon_interp;
            amrex::Real const iy2 = has_sub ? eiy(i,j,k,2) : amrex::Real(0.0);
            Ey(i, j, k) = alpha * Ey(i, j, k)
                        + beta0 * ( iy0*Dx0 + iy1*Dy0 + iy2*Dz0 )
                        - beta0 * iy1 * jy(i, j, k);
        });

        // Ez update
        amrex::ParallelFor(tez, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // Skip field push in the embedded boundaries
            if (update_Ez_arr && update_Ez_arr(i, j, k) == 0) { return; }
            // sigma > 0 cells are handled by the exponential step below; leaving them
            // untouched here keeps E^n available for the matrix product.
            // node sigma drives both the skip and the damping (§4.1)
            amrex::Real const sigma_node = has_sig_mat ? sn2(i, j, k, 0) : amrex::Real(0.0);
            if (has_sig_mat && sigma_node > amrex::Real(0.0)) { return; }   // exponential step handles it

            // Interpolate conductivity, sigma, to Ez position on the grid
            // sigma-free cells must not inherit damping from a neighbouring conductor
            // through the face interpolation: that made the vacuum side of a sigma
            // interface decay and locked E (and B) to zero inside the conductor
            // (c03-E6, 技术路径 §7.11). With sigma = 0 everywhere this is a no-op.
            amrex::Real const sigma_interp = sigma_node;   // node sigma, not a cell interpolation
            // Interpolated permittivity, epsilon, to Ez position on the grid
            amrex::Real const epsilon_interp = ablastr::coarsen::sample::Interp(eps_arr, epsilon_stag,
                                                                                Ez_stag, macro_cr, i, j, k, scomp);
            const amrex::Real alpha = T_MacroAlgo::alpha( sigma_interp, epsilon_interp, dt);
            const amrex::Real beta = T_MacroAlgo::beta( sigma_interp, epsilon_interp, dt);

            amrex::Real const Dz0 =   T_Algo::DownwardDx(Hy, coefs_x, n_coefs_x, i, j, k,0)
                                    - T_Algo::DownwardDy(Hx, coefs_y, n_coefs_y, i, j, k,0);
            amrex::Real const Dx0 = - T_Algo::DownwardDz(Hy, coefs_z, n_coefs_z, i, j, k,0)
                                    + T_Algo::DownwardDy(Hz, coefs_y, n_coefs_y, i, j, k,0);
            amrex::Real const Dy0 =   T_Algo::DownwardDz(Hx, coefs_z, n_coefs_z, i, j, k,0)
                                    - T_Algo::DownwardDx(Hz, coefs_x, n_coefs_x, i, j, k,0);
            amrex::Real const beta0 = beta * epsilon_interp;
            amrex::Real const iz0 = has_sub ? eiz(i,j,k,0) : amrex::Real(0.0);
            amrex::Real const iz1 = has_sub ? eiz(i,j,k,1) : amrex::Real(0.0);
            amrex::Real const iz2 = has_sub ? eiz(i,j,k,2) : amrex::Real(1.0)/epsilon_interp;
            Ez(i, j, k) = alpha * Ez(i, j, k)
                        + beta0 * ( iz0*Dx0 + iz1*Dy0 + iz2*Dz0 )
                        - beta0 * iz2 * jz(i, j, k);
        });
        // P2: exact exponential step for the conductive cells:
        //     E^{n+1} = M1 E^n + M2 (curlH - J)
        // Same curl expressions and the same sign for J as the Yee branch above.

    }

    // Exponential step on the NODE layout: one top-level loop per E component (AMReX
    // forbids nested MFIters, so this cannot live inside the BE loop above). Each node
    // is written in its own component layout, so BOTH grid copies of a shared node get
    // the same value -- this is what makes the result independent of where the sigma
    // interface sits relative to the grid boundaries (V1, see 技术路径 §7.11 / §4.1).
    if (macroscopic_properties->has_sigma_matrices()) {
        Real const * const cxs = m_stencil_coefs_x.dataPtr();
        auto const ncx = static_cast<int>(m_stencil_coefs_x.size());
        Real const * const cys = m_stencil_coefs_y.dataPtr();
        auto const ncy = static_cast<int>(m_stencil_coefs_y.size());
        Real const * const czs = m_stencil_coefs_z.dataPtr();
        auto const ncz = static_cast<int>(m_stencil_coefs_z.size());
        for (int c = 0; c < 3; ++c) {
            auto const& m1n = macroscopic_properties->get_sigma_m1_node_mf(c);
            auto const& m2n = macroscopic_properties->get_sigma_m2_node_mf(c);
            auto const& snn = macroscopic_properties->get_sigma_node_mf(c);
            for (MFIter mfi_n(*Efield[c], TilingIfNotGPU()); mfi_n.isValid(); ++mfi_n) {
                Box const& tb = mfi_n.validbox();
                auto const mun = macroscopic_properties->getmu_mf().array(mfi_n);
                auto const b0 = Bfield[0]->array(mfi_n);
                auto const b1 = Bfield[1]->array(mfi_n);
                auto const b2 = Bfield[2]->array(mfi_n);
                auto const m1 = m1n.const_array(mfi_n);
                auto const m2 = m2n.const_array(mfi_n);
                auto const sn = snn.const_array(mfi_n);
                auto const ex = Efield[0]->const_array(mfi_n);
                auto const ey = Efield[1]->const_array(mfi_n);
                auto const ez = Efield[2]->const_array(mfi_n);
                auto const jx = Jfield[0]->const_array(mfi_n);
                auto const jy = Jfield[1]->const_array(mfi_n);
                auto const jz = Jfield[2]->const_array(mfi_n);
                auto Ev = Efield[c]->array(mfi_n);
                amrex::Array4<int> upd_c{};
                if (eb_update_E[c] != nullptr) { upd_c = eb_update_E[c]->array(mfi_n); }
                FieldAccessorMacroscopic const Hx(b0, mun);
                FieldAccessorMacroscopic const Hy(b1, mun);
                FieldAccessorMacroscopic const Hz(b2, mun);
                amrex::ParallelFor(tb, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (sn(i, j, k, 0) <= amrex::Real(0.0)) { return; }
                    amrex::Real const cx = - T_Algo::DownwardDz(Hy, czs, ncz, i, j, k,0)
                                           + T_Algo::DownwardDy(Hz, cys, ncy, i, j, k,0);
                    amrex::Real const cy =   T_Algo::DownwardDz(Hx, czs, ncz, i, j, k,0)
                                           - T_Algo::DownwardDx(Hz, cxs, ncx, i, j, k,0);
                    amrex::Real const cz =   T_Algo::DownwardDx(Hy, cxs, ncx, i, j, k,0)
                                           - T_Algo::DownwardDy(Hx, cys, ncy, i, j, k,0);
                    amrex::Real const rhs[3] = { cx - jx(i, j, k), cy - jy(i, j, k), cz - jz(i, j, k) };
                    amrex::Real const eold[3] = { ex(i, j, k), ey(i, j, k), ez(i, j, k) };
                    amrex::Real en = amrex::Real(0.0);
                    for (int d = 0; d < 3; ++d) {
                        en += m1(i, j, k, 3*c + d) * eold[d]
                            + m2(i, j, k, 3*c + d) * rhs[d];
                    }
                    if (!(upd_c && upd_c(i, j, k) == 0)) { Ev(i, j, k) = en; }
                });
            }
        }
    }
}

#endif // corresponds to ifndef WARPX_DIM_RZ
