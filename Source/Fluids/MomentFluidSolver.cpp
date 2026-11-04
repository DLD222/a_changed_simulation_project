/* Copyright (c) 2026
 *
 * See MomentFluidSolver.H for the contract.  Prototype: 2D XZ, single box,
 * periodic, host-side.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Fluids/MomentFluidSolver.H"

#include "MomentsSolver/FV/MomentFVNDDriver.H"
#include "MomentsSolver/Coupling/YeeDeposition.H"
#include "MomentsSolver/Equations/Euler5Eqn.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Loop.H>
#include <AMReX_BLassert.H>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace amrex;

namespace {
using MEq = moments::Euler5;
constexpr int MNEQ = MEq::num_eqs;   // 5
constexpr Real s_pi = Real(3.14159265358979323846264338327950288);

/** Read one real from `fluids.moment.<name>` with a default. */
Real get_real (char const * name, Real dflt)
{
    Real v = dflt;
    ParmParse const pp("fluids.moment");
    pp.query(name, v);
    return v;
}
} // namespace

bool MomentFluidSolver::IsEnabled ()
{
    static bool const enabled = [] () {
        int flag = 0;
        ParmParse const pp("fluids.moment");
        pp.query("enabled", flag);
        return flag != 0;
    }();
    return enabled;
}

void MomentFluidSolver::InitData (Geometry const & geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == 2,
        "fluids.moment: only 2D XZ is supported by this prototype");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.isPeriodic(0) && geom.isPeriodic(1),
        "fluids.moment: both directions must be periodic in this prototype");

    m_q     = get_real("q",         m_q);
    m_m     = get_real("m",         m_m);
    m_gamma = get_real("gamma",     m_gamma);
    m_eps0  = get_real("eps0",      m_eps0);
    m_amp   = get_real("amplitude", m_amp);
    m_p0    = get_real("pressure",  m_p0);

    BoxArray const ba(geom.Domain());
    DistributionMapping const dm(ba);
    m_U = std::make_unique<MultiFab>(ba, dm, MNEQ, 0);

    Real const dx  = geom.CellSize(0);
    Real const lo0 = geom.ProbLo(0);
    Real const k   = Real(2.0) * s_pi / geom.ProbLength(0);
    Real const gam = m_gamma;

    for (MFIter mfi(*m_U); mfi.isValid(); ++mfi) {
        Box const & bx = mfi.validbox();
        auto const & u = m_U->array(mfi);
        Real const amp = m_amp;
        Real const p0  = m_p0;
        amrex::LoopOnCpu(bx, [&] (int i, int j, int kk) {
            (void)kk;
            Real const x   = lo0 + (i + Real(0.5)) * dx;
            Real const rho = Real(1.0) + amp * std::cos(k * x);
            u(i,j,kk,0) = rho;
            u(i,j,kk,1) = Real(0.0);
            u(i,j,kk,2) = Real(0.0);
            u(i,j,kk,3) = Real(0.0);
            u(i,j,kk,4) = p0 / (gam - Real(1.0));
        });
    }
    m_U->FillBoundary(geom.periodicity());
}

void MomentFluidSolver::AdvanceAndDeposit (Real dt, MultiFab* J[3],
                                           Geometry const & geom, int lev)
{
    if (lev != 0) { return; }
    if (m_U == nullptr) { InitData(geom); }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_U->boxArray().size() == 1,
        "fluids.moment: this prototype needs a single box (raise amr.max_grid_size)");

    Box const & bx = m_U->boxArray()[0];
    int const nx = static_cast<int>(bx.length(0));
    int const nz = static_cast<int>(bx.length(1));
    Real const dx = geom.CellSize(0);
    Real const dz = geom.CellSize(1);
    int const i0 = bx.smallEnd(0);
    int const j0 = bx.smallEnd(1);

    // ---- gather into a contiguous host buffer ----
    std::vector<Real> U(static_cast<std::size_t>(nx) * nz * MNEQ);
    {
        auto const & u = m_U->const_array(0);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            (void)k;
            std::size_t const c = static_cast<std::size_t>(i - i0) * nz + (j - j0);
            for (int v = 0; v < MNEQ; ++v) { U[c * MNEQ + v] = u(i,j,k,v); }
        });
    }

    // ---- advance, recording the RK-weighted mass flux ----
    moments::FVND<MEq, 2> fv(nx, nz, 1, dx, dz, Real(1.0));
    MEq::Context ctx;
    ctx.gamma = m_gamma;
    moments::BC mbc[3] = {moments::BC::Periodic, moments::BC::Periodic,
                          moments::BC::Periodic};
    moments::YeeDeposit<2> dep;
    moments::deposit_conservative<MEq, 2>(fv, U.data(), dt, ctx, mbc,
                                          m_q, m_m, dep);

    // ---- scatter back ----
    {
        auto const & u = m_U->array(0);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            (void)k;
            std::size_t const c = static_cast<std::size_t>(i - i0) * nz + (j - j0);
            for (int v = 0; v < MNEQ; ++v) { u(i,j,k,v) = U[c * MNEQ + v]; }
        });
    }

    // ---- add the current to the Yee faces ----
    // dir 0 -> x faces (J[0]); dir 1 -> z faces (J[2]).  J[1] (y) is untouched.
    Real const scale = m_q / m_m;
    for (int dir = 0; dir < 2; ++dir) {
        MultiFab * jmf = (dir == 0) ? J[0] : J[2];
        if (jmf == nullptr) { continue; }
        for (MFIter mfi(*jmf); mfi.isValid(); ++mfi) {
            Box const & jbx = mfi.validbox();
            auto const & j = jmf->array(mfi);
            amrex::LoopOnCpu(jbx, [&] (int i, int j2, int k) {
                (void)k;
                std::size_t const c = static_cast<std::size_t>(i - i0) * nz + (j2 - j0);
                j(i,j2,k) += scale * dep.Jlo[dir][c];
            });
        }
        jmf->FillBoundary(geom.periodicity());
    }
}

Real MomentFluidSolver::GaussResidual (MultiFab const & Ex,
                                       MultiFab const & Ez,
                                       Geometry const & geom) const
{
    if (m_U == nullptr) { return Real(0.0); }

    Real const dx = geom.CellSize(0);
    Real const dz = geom.CellSize(1);
    Real const scale = m_q / m_m;
    Real worst = Real(0.0);

    for (MFIter mfi(*m_U); mfi.isValid(); ++mfi) {
        Box const & bx = mfi.validbox();
        auto const & u  = m_U->const_array(mfi);
        auto const & ex = Ex.const_array(mfi);
        auto const & ez = Ez.const_array(mfi);
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                Real const divE = (ex(i+1,j,0) - ex(i,j,0)) / dx
                                + (ez(i,j+1,0) - ez(i,j,0)) / dz;
                // the fluid charge plus the neutralising background (+1)
                Real const rho_q = scale * u(i,j,0,0) + Real(1.0);
                worst = std::max(worst, std::abs(divE - rho_q / m_eps0));
            }
        }
    }
    return worst;
}
