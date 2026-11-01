#include <fstream>
#include "MacroscopicProperties.H"
#include "SubpixelTensor.H"
#include "WarpX.H"
#include <AMReX_EBFabFactory.H>

#include "Fields.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_Array4.H>
#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Parser.H>

#include <AMReX_BaseFwd.H>

#include <memory>
#include <sstream>

using namespace amrex;

MacroscopicProperties::MacroscopicProperties ()
{
    ReadParameters();
}

void
MacroscopicProperties::ReadParameters ()
{
    const ParmParse pp_macroscopic("macroscopic");
    // Since macroscopic maxwell solve is turned on,
    // user-defined sigma, mu, and epsilon are queried.
    // The vacuum values are used as default for the macroscopic parameters
    // with a warning message to the user to indicate that no value was specified.

    // Query input for material conductivity, sigma.
    bool sigma_specified = false;
    if (utils::parser::queryWithParser(pp_macroscopic, "sigma", m_sigma)) {
        m_sigma_s = "constant";
        sigma_specified = true;
    }
    if (pp_macroscopic.query("sigma_function(x,y,z)", m_str_sigma_function) ) {
        m_sigma_s = "parse_sigma_function";
        sigma_specified = true;
    }
    if (!sigma_specified) {
        std::stringstream warnMsg;
        warnMsg << "Material conductivity is not specified. Using default vacuum value of " <<
            m_sigma << " in the simulation.";
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            warnMsg.str());
    }
    // initialization of sigma (conductivity) with parser
    if (m_sigma_s == "parse_sigma_function") {
        utils::parser::Store_parserString(
            pp_macroscopic, "sigma_function(x,y,z)", m_str_sigma_function);
        m_sigma_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(m_str_sigma_function,{"x","y","z"}));
    }

    // Query input for material permittivity, epsilon.
    bool epsilon_specified = false;
    if (utils::parser::queryWithParser(pp_macroscopic, "epsilon", m_epsilon)) {
        m_epsilon_s = "constant";
        epsilon_specified = true;
    }
    if (pp_macroscopic.query("epsilon_function(x,y,z)", m_str_epsilon_function) ) {
        m_epsilon_s = "parse_epsilon_function";
        epsilon_specified = true;
    }
    // A tensor key also counts as "specified": the scalar `epsilon` is then ignored.
    // Detecting it here (before the tensor block below) suppresses the misleading
    // "using default vacuum value" warning in tensor mode.
    bool tensor_key = false;
    {
        char const* tk[12] = {"epsilon_tensor_xx","epsilon_tensor_yy","epsilon_tensor_zz",
                              "epsilon_tensor_xy","epsilon_tensor_xz","epsilon_tensor_yz",
                              "epsilon2_tensor_xx","epsilon2_tensor_yy","epsilon2_tensor_zz",
                              "epsilon2_tensor_xy","epsilon2_tensor_xz","epsilon2_tensor_yz"};
        for (int i = 0; i < 12; ++i) { tensor_key = tensor_key || pp_macroscopic.contains(tk[i]); }
    }
    if (!epsilon_specified && !tensor_key) {
        std::stringstream warnMsg;
        warnMsg << "Material permittivity is not specified. Using default vacuum value of " <<
            m_epsilon << " in the simulation.";
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            warnMsg.str());
    }

    // Anisotropic permittivity: 6 independent parser entries per material
    // (xx,yy,zz,xy,xz,yz), absolute SI, same units as `macroscopic.epsilon`.
    // Separate keys (instead of one comma/space separated string) because parser
    // expressions may themselves contain commas or spaces.
    {
        char const* k1[6] = {"epsilon_tensor_xx","epsilon_tensor_yy","epsilon_tensor_zz",
                             "epsilon_tensor_xy","epsilon_tensor_xz","epsilon_tensor_yz"};
        char const* k2[6] = {"epsilon2_tensor_xx","epsilon2_tensor_yy","epsilon2_tensor_zz",
                             "epsilon2_tensor_xy","epsilon2_tensor_xz","epsilon2_tensor_yz"};
        // Presence is decided by the KEY, not by a successful parse: a key that exists
        // but fails to parse must abort, never silently fall back to the scalar path.
        // (This is the failure mode that produced Ez=0 / v=0.808c with a "normal" run.)
        bool key1 = false, key2 = false;
        for (int i = 0; i < 6; ++i) {
            key1 = key1 || pp_macroscopic.contains(k1[i]);
            key2 = key2 || pp_macroscopic.contains(k2[i]);
        }
        amrex::Real v1[6] = {0,0,0,0,0,0};
        amrex::Real v2[6] = {0,0,0,0,0,0};
        int n1 = 0, n2 = 0;
        for (int i = 0; i < 6; ++i) {
            // Per-component parser: a constant ("1.0") or a spatial function ("1+x").
            // makeParser handles both; the value at (0,0,0) is kept for the asserts/print.
            if (pp_macroscopic.query(k1[i], m_str_eps_tensor[i])) {
                m_eps_tensor_parser[i] = std::make_unique<amrex::Parser>(
                    utils::parser::makeParser(m_str_eps_tensor[i], {"x","y","z"}));
                v1[i] = m_eps_tensor_parser[i]->compile<3>()(amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0));
                ++n1;
            }
            if (pp_macroscopic.query(k2[i], m_str_eps2_tensor[i])) {
                m_eps2_tensor_parser[i] = std::make_unique<amrex::Parser>(
                    utils::parser::makeParser(m_str_eps2_tensor[i], {"x","y","z"}));
                v2[i] = m_eps2_tensor_parser[i]->compile<3>()(amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0));
                ++n2;
            }
        }
        if (key1) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(n1 > 0,
                "macroscopic: epsilon_tensor_* is present but no component parsed; "
                "refusing to fall back to the scalar epsilon");
            bool nz = false;
            for (int i = 0; i < 6; ++i) { nz = nz || (v1[i] != amrex::Real(0.0)); }
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(nz,
                "macroscopic: epsilon_tensor_* parsed to all zeros");
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!epsilon_specified,
                "macroscopic: give either `epsilon` or `epsilon_tensor_*`, not both");
        }
        if (key2) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(n2 > 0,
                "macroscopic: epsilon2_tensor_* is present but no component parsed");
        }
        for (int i = 0; i < 6; ++i) {
            m_eps_tensor[i]  = v1[i];
            m_eps2_tensor[i] = v2[i];
        }
        m_has_eps_tensor  = key1;
        m_has_eps2_tensor = key2;
        if (key1) {
            amrex::Print() << "[MACRO] eps1 tensor: xx=" << v1[0] << " yy=" << v1[1]
                           << " zz=" << v1[2] << " xy=" << v1[3] << " xz=" << v1[4]
                           << " yz=" << v1[5] << "\n";
        }
        // Material 1 may be scalar or tensor, independently of material 2. A scalar
        // eps1 + tensor eps2 (or the reverse) is a legitimate interface configuration;
        // the Kottke rows and fill_sigma_mats each fall back to diag(eps) per material.
    }


    // Dielectric interface (independent of the PEC embedded boundaries).
    // phi > 0 marks material 1; material 2 is described by epsilon2.
    pp_macroscopic.query("epsilon2", m_epsilon2);
    if (pp_macroscopic.query("interface_function(x,y,z)", m_str_interface_function)) {
        utils::parser::Store_parserString(
            pp_macroscopic, "interface_function(x,y,z)", m_str_interface_function);
        m_interface_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(m_str_interface_function,{"x","y","z"}));
        m_has_interface = true;
    }

    // initialization of epsilon (permittivity) with parser
    if (m_epsilon_s == "parse_epsilon_function") {
        utils::parser::Store_parserString(
            pp_macroscopic, "epsilon_function(x,y,z)", m_str_epsilon_function);
        m_epsilon_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(m_str_epsilon_function,{"x","y","z"}));
    }

    // Query input for material permeability, mu.
    bool mu_specified = false;
    if (utils::parser::queryWithParser(pp_macroscopic, "mu", m_mu)) {
        m_mu_s = "constant";
        mu_specified = true;
    }
    if (pp_macroscopic.query("mu_function(x,y,z)", m_str_mu_function) ) {
        m_mu_s = "parse_mu_function";
        mu_specified = true;
    }
    if (!mu_specified) {
        std::stringstream warnMsg;
        warnMsg << "Material permeability is not specified. Using default vacuum value of " <<
            m_mu << " in the simulation.";
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            warnMsg.str());
    }

    // initialization of mu (permeability) with parser
    if (m_mu_s == "parse_mu_function") {
        utils::parser::Store_parserString(
            pp_macroscopic, "mu_function(x,y,z)", m_str_mu_function);
        m_mu_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(m_str_mu_function,{"x","y","z"}));
    }

}

void
MacroscopicProperties::AllocateLevelMFs (
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dmap,
    const amrex::IntVect& ng_EB_alloc )
{
    // Need >=1 ghost cell so that i-1 / j-1 reads at the domain boundary are in-bounds
    // (seed/update/diagnose of rho_medium and the per-component Interp read neighbours).
    amrex::IntVect ng_props = ng_EB_alloc;
    ng_props = amrex::max(ng_props, amrex::IntVect::TheUnitVector());
    // sigma is cell-centered MultiFab
    m_sigma_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 1, ng_props);
    // epsilon is cell-centered MultiFab
    m_eps_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 1, ng_props);
    // mu is cell-centered MultiFab
    m_mu_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 1, ng_props);
    // Cell-centred anisotropic permittivity tensors (6 comps each), present only when
    // epsilon_tensor_* / epsilon2_tensor_* are given.
    if (m_has_eps_tensor) {
        m_eps_tensor_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 6, ng_props);
        m_eps_tensor_mf->setVal(amrex::Real(0.0));
    }
    if (m_has_eps2_tensor) {
        m_eps2_tensor_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 6, ng_props);
        m_eps2_tensor_mf->setVal(amrex::Real(0.0));
    }
}

void
MacroscopicProperties::InitData (
    const amrex::Geometry& geom,
    const amrex::IntVect& Ex_stag,
    const amrex::IntVect& Ey_stag,
    const amrex::IntVect& Ez_stag)
{
    // Guard against the silent downgrade that bit c03-E3: with an interface_function
    // but subpixel off, the subpixel rows are never built and epsilon2 / tensor
    // material 2 are simply ignored (the run then looks exactly like the vacuum run).
    if (m_has_interface && !WarpX::macroscopic_subpixel) {
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            "interface_function is set but algo.macroscopic_subpixel=0: the subpixel "
            "rows are NOT built and epsilon2 / tensor material 2 will be ignored.");
    }
    // Initialize sigma
    if (m_sigma_s == "constant") {

        m_sigma_mf->setVal(m_sigma);

    } else if (m_sigma_s == "parse_sigma_function") {

        InitializeMacroMultiFabUsingParser(m_sigma_mf.get(), m_sigma_parser->compile<3>(),
            geom.CellSizeArray(), geom.ProbDomain());
    }
    // Initialize epsilon
    if (m_epsilon_s == "constant") {

        m_eps_mf->setVal(m_epsilon);

    } else if (m_epsilon_s == "parse_epsilon_function") {

        InitializeMacroMultiFabUsingParser(m_eps_mf.get(), m_epsilon_parser->compile<3>(),
        geom.CellSizeArray(), geom.ProbDomain());

    }
    // In the Maxwell solver, `epsilon` is used in the denominator.
    // Therefore, it needs to be strictly positive
    bool const local=true;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( m_eps_mf->min(0,0,local) > 0,
    "WarpX encountered zero or negative values for the relative permittivity `epsilon`. Please check the initialization of `epsilon`.");

    // Initialize mu
    if (m_mu_s == "constant") {

        m_mu_mf->setVal(m_mu);

    } else if (m_mu_s == "parse_mu_function") {

        InitializeMacroMultiFabUsingParser(m_mu_mf.get(), m_mu_parser->compile<3>(),
            geom.CellSizeArray(), geom.ProbDomain());

    }

    // Fill cell-centred anisotropic permittivity tensors (6 comps each) from the
    // per-component parsers (constant or (x,y,z) function). A component without a
    // parser stays 0, matching the previous GpuArray default.
    if (m_has_eps_tensor || m_has_eps2_tensor) {
        auto const dxa = geom.CellSizeArray();
        amrex::Real const lox = geom.ProbDomain().lo(0);
        amrex::Real const loz = geom.ProbDomain().lo(AMREX_SPACEDIM-1);
        auto fill_one = [&](amrex::MultiFab& mf, int c, amrex::Parser const* p) {
            if (p == nullptr) { return; }
            amrex::ParserExecutor<3> const pex = p->compile<3>();
            for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box const& tb = mfi.tilebox(mf.ixType().toIntVect(), mf.nGrowVect());
                auto const& fab = mf.array(mfi);
                amrex::ParallelFor(tb, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    amrex::Real const x = (amrex::Real(i) + amrex::Real(0.5)) * dxa[0] + lox;
                    amrex::Real const z = (amrex::Real(j) + amrex::Real(0.5)) * dxa[1] + loz;
                    fab(i,j,k,c) = pex(x, amrex::Real(0.0), z);
                });
            }
        };
        for (int c = 0; c < 6; ++c) {
            if (m_has_eps_tensor) { fill_one(*m_eps_tensor_mf, c, m_eps_tensor_parser[c].get()); }
            if (m_has_eps2_tensor) { fill_one(*m_eps2_tensor_mf, c, m_eps2_tensor_parser[c].get()); }
        }
    }

    // Guard: sigma > 0 with subpixel off silently falls back to the first-order
    // backward-Euler integrator. The exponential M1/M2 are only built inside
    // InitSubpixelRows, which is never called when macroscopic_subpixel = 0;
    // at a*dt ~ O(0.1..1) the two integrators differ visibly (see c03-E6).
    if (!WarpX::macroscopic_subpixel && m_sigma_mf && m_sigma_mf->max(0) > amrex::Real(0.0)) {
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            "macroscopic.sigma > 0 but algo.macroscopic_subpixel=0: the E update "
            "falls back to the first-order backward-Euler integrator instead of "
            "the exact exponential step.");
    }

    // Fill the ghost cells once: the parser/constant init already wrote the analytic
    // continuation into the ghost region, but FillBoundary exchanges the periodic
    // boundary values so that i-1 / j-1 reads at the domain edge are the wrapped values.
    // eps/sigma/mu are static, so this is a one-time cost, not per-step.
    m_eps_mf->FillBoundary(geom.periodicity());
    m_sigma_mf->FillBoundary(geom.periodicity());
    m_mu_mf->FillBoundary(geom.periodicity());

    const amrex::IntVect sigma_stag = m_sigma_mf->ixType().toIntVect();
    const amrex::IntVect epsilon_stag = m_eps_mf->ixType().toIntVect();
    const amrex::IntVect mu_stag = m_mu_mf->ixType().toIntVect();


    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        sigma_IndexType[idim]   = sigma_stag[idim];
        epsilon_IndexType[idim] = epsilon_stag[idim];
        mu_IndexType[idim]      = mu_stag[idim];
        Ex_IndexType[idim]      = Ex_stag[idim];
        Ey_IndexType[idim]      = Ey_stag[idim];
        Ez_IndexType[idim]      = Ez_stag[idim];
        macro_cr_ratio[idim]    = 1;
    }
// Initialize components that were left uninitialized by the above loop
#if defined(WARPX_DIM_1D_Z)
        sigma_IndexType[1]   = 0;  sigma_IndexType[2]   = 0;
        epsilon_IndexType[1] = 0;  epsilon_IndexType[2] = 0;
        mu_IndexType[1]      = 0;  mu_IndexType[2]      = 0;
        Ex_IndexType[1]      = 0;  Ex_IndexType[2]      = 0;
        Ey_IndexType[1]      = 0;  Ey_IndexType[2]      = 0;
        Ez_IndexType[1]      = 0;  Ez_IndexType[2]      = 0;
        macro_cr_ratio[1]    = 1;  macro_cr_ratio[2]    = 1;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        sigma_IndexType[2]   = 0;
        epsilon_IndexType[2] = 0;
        mu_IndexType[2]      = 0;
        Ex_IndexType[2]      = 0;
        Ey_IndexType[2]      = 0;
        Ez_IndexType[2]      = 0;
        macro_cr_ratio[2]    = 1;
#endif
}

void
MacroscopicProperties::InitializeMacroMultiFabUsingParser (
                       amrex::MultiFab *macro_mf,
                       amrex::ParserExecutor<3> const& macro_parser,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_lev,
                       const amrex::RealBox& prob_domain_lev)
{
    const amrex::IntVect iv = macro_mf->ixType().toIntVect();
    for ( amrex::MFIter mfi(*macro_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // Initialize ghost cells in addition to valid cells

        const amrex::Box& tb = mfi.tilebox( iv, macro_mf->nGrowVect());
        amrex::Array4<amrex::Real> const& macro_fab =  macro_mf->array(mfi);
        amrex::ParallelFor (tb,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Shift x, y, z position based on index type

#if defined(WARPX_DIM_1D_Z)
                const amrex::Real x = 0._rt;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - iv[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real z = i * dx_lev[0] + prob_domain_lev.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const amrex::Real fac_x = (1._rt - iv[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i * dx_lev[0] + prob_domain_lev.lo(0) + fac_x;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real z = j * dx_lev[1] + prob_domain_lev.lo(1) + fac_z;
#else
                const amrex::Real fac_x = (1._rt - iv[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i * dx_lev[0] + prob_domain_lev.lo(0) + fac_x;
                const amrex::Real fac_y = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real y = j * dx_lev[1] + prob_domain_lev.lo(1) + fac_y;
                const amrex::Real fac_z = (1._rt - iv[2]) * dx_lev[2] * 0.5_rt;
                const amrex::Real z = k * dx_lev[2] + prob_domain_lev.lo(2) + fac_z;
#endif
                // initialize the macroparameter
                macro_fab(i,j,k) = macro_parser(x,y,z);
        });

    }
}

void
MacroscopicProperties::InitSubpixelRows (
    const amrex::Geometry& geom,
    const amrex::EBFArrayBoxFactory* eb_fact,
    const amrex::IntVect& Ex_stag,
    const amrex::IntVect& Ey_stag,
    const amrex::IntVect& Ez_stag,
    amrex::Real dt)
{
    amrex::ignore_unused(eb_fact);
    m_sigma_dt = dt;

    auto const& ba = m_eps_mf->boxArray();
    auto const& dm = m_eps_mf->DistributionMap();

    // P2: precompute M1 = exp(-A dt), M2 = (1-exp(-A dt)) A^-1 eps^-1 for sigma>0
    // cells only (see 技术路径 §7.7). Called from every branch that fills the rows.
    auto fill_sigma_mats = [&]() {
        amrex::Real const smax = m_sigma_mf ? m_sigma_mf->max(0) : amrex::Real(0.0);
        if (smax <= amrex::Real(0.0)) { return; }
        m_sigma_m1 = std::make_unique<amrex::MultiFab>(ba, dm, 9, 0);
        m_sigma_m2 = std::make_unique<amrex::MultiFab>(ba, dm, 9, 0);
        m_sigma_m1->setVal(amrex::Real(0.0));
        m_sigma_m2->setVal(amrex::Real(0.0));
        amrex::MultiFab const& sigm = *m_sigma_mf;
        amrex::MultiFab const& epsm = *m_eps_mf;
        // eps here must be the eps of THIS cell, not material 1: m_eps_mf stores only
        // material 1, so using it everywhere made A = eps^-1 sigma wrong by O(1) in
        // material 2 (4x when eps2/eps1 = 4). Reuse the same phi sub-sampling that the
        // Kottke rows use, so both agree on which material each cell is.
        amrex::Real const eps2s = (m_epsilon2 > amrex::Real(0.0)) ? m_epsilon2 : m_epsilon;
        bool const iface_on = (m_has_interface && (m_interface_parser != nullptr));
        amrex::Real const plo_x = geom.ProbDomain().lo(0);
        amrex::Real const plo_z = geom.ProbDomain().lo(AMREX_SPACEDIM-1);
        amrex::Real const dxs = geom.CellSize(0);
        amrex::Real const dzs = geom.CellSize(AMREX_SPACEDIM-1);
        auto const phi_s = iface_on ? m_interface_parser->compile<3>() : amrex::ParserExecutor<3>{};
        constexpr int nsub2 = 4;
        amrex::Real const inv_ns = amrex::Real(1.0) / static_cast<amrex::Real>(nsub2);
        bool const has_tensor_e  = m_has_eps_tensor;
        bool const has_tensor_e2 = m_has_eps2_tensor;
        for (amrex::MFIter mfi(*m_sigma_m1, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.validbox();
            auto const& m1a = m_sigma_m1->array(mfi);
            auto const& m2a = m_sigma_m2->array(mfi);
            auto const& sa  = sigm.const_array(mfi);
            auto const& ea  = epsm.const_array(mfi);
            auto const& et1 = has_tensor_e  ? m_eps_tensor_mf->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            auto const& et2 = has_tensor_e2 ? m_eps2_tensor_mf->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::Real const sv = sa(i,j,k);
                if (sv <= amrex::Real(0.0)) { return; }
                // fill fraction of material 1 (same phi sub-sampling as the Kottke rows)
                amrex::Real f1 = amrex::Real(1.0);
                if (iface_on) {
                    amrex::Real const xc = plo_x + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dxs;
                    amrex::Real const zc = plo_z + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dzs;
                    int n1 = 0;
                    for (int a = 0; a < nsub2; ++a) {
                        for (int b = 0; b < nsub2; ++b) {
                            amrex::Real const x = xc + (static_cast<amrex::Real>(a) + amrex::Real(0.5)) * inv_ns * dxs - amrex::Real(0.5) * dxs;
                            amrex::Real const z = zc + (static_cast<amrex::Real>(b) + amrex::Real(0.5)) * inv_ns * dzs - amrex::Real(0.5) * dzs;
                            if (phi_s(x, amrex::Real(0.0), z) > amrex::Real(0.0)) { ++n1; }
                        }
                    }
                    f1 = static_cast<amrex::Real>(n1) * inv_ns * inv_ns;
                }
                // A = eps^-1 sigma: eps is the local tensor (6 comps: xx,yy,zz,xy,xz,yz),
                // f-weighted at an interface; material 2 falls back to scalar eps2 on the
                // diagonal when only eps1 is a tensor.
                warpx::subpixel::Mat3 E{};
                amrex::Real const w2 = amrex::Real(1.0) - f1;
                if (has_tensor_e || has_tensor_e2) {
                    // either material may be scalar (diag) or tensor (6 comps)
                    amrex::Real const e1_xx = has_tensor_e  ? et1(i,j,k,0) : ea(i,j,k);
                    amrex::Real const e1_yy = has_tensor_e  ? et1(i,j,k,1) : ea(i,j,k);
                    amrex::Real const e1_zz = has_tensor_e  ? et1(i,j,k,2) : ea(i,j,k);
                    amrex::Real const e1_xy = has_tensor_e  ? et1(i,j,k,3) : amrex::Real(0.0);
                    amrex::Real const e1_xz = has_tensor_e  ? et1(i,j,k,4) : amrex::Real(0.0);
                    amrex::Real const e1_yz = has_tensor_e  ? et1(i,j,k,5) : amrex::Real(0.0);
                    amrex::Real const e2_xx = has_tensor_e2 ? et2(i,j,k,0) : eps2s;
                    amrex::Real const e2_yy = has_tensor_e2 ? et2(i,j,k,1) : eps2s;
                    amrex::Real const e2_zz = has_tensor_e2 ? et2(i,j,k,2) : eps2s;
                    amrex::Real const e2_xy = has_tensor_e2 ? et2(i,j,k,3) : amrex::Real(0.0);
                    amrex::Real const e2_xz = has_tensor_e2 ? et2(i,j,k,4) : amrex::Real(0.0);
                    amrex::Real const e2_yz = has_tensor_e2 ? et2(i,j,k,5) : amrex::Real(0.0);
                    E(0,0) = f1*e1_xx + w2*e2_xx;
                    E(1,1) = f1*e1_yy + w2*e2_yy;
                    E(2,2) = f1*e1_zz + w2*e2_zz;
                    E(0,1) = E(1,0) = f1*e1_xy + w2*e2_xy;
                    E(0,2) = E(2,0) = f1*e1_xz + w2*e2_xz;
                    E(1,2) = E(2,1) = f1*e1_yz + w2*e2_yz;
                } else {
                    amrex::Real const ev = f1 * ea(i,j,k) + w2 * eps2s;
                    E(0,0) = ev; E(1,1) = ev; E(2,2) = ev;
                }
                warpx::subpixel::Mat3 S{};
                S(0,0) = sv; S(1,1) = sv; S(2,2) = sv;
                warpx::subpixel::Mat3 MM1{}, MM2{};
                warpx::subpixel::sigma_exponential_matrices(E, S, dt, MM1, MM2);
                for (int c = 0; c < 9; ++c) { m1a(i,j,k,c) = MM1.m[c]; m2a(i,j,k,c) = MM2.m[c]; }
            });
        }
        m_sigma_mat_ready = true;
        amrex::Print() << "[SIGMA-EXP] sigma_max=" << smax << " dt=" << dt
                       << "  M1/M2 filled where sigma>0\n";
    };
    fill_sigma_mats();   // single exit point: independent of the rows branch
    int const ng = m_eps_mf->nGrow();

    m_eps_inv_x = std::make_unique<amrex::MultiFab>(amrex::convert(ba, Ex_stag), dm, 3, ng);
    m_eps_inv_y = std::make_unique<amrex::MultiFab>(amrex::convert(ba, Ey_stag), dm, 3, ng);
    m_eps_inv_z = std::make_unique<amrex::MultiFab>(amrex::convert(ba, Ez_stag), dm, 3, ng);
    m_eps_inv_x->setVal(amrex::Real(0.0));
    m_eps_inv_y->setVal(amrex::Real(0.0));
    m_eps_inv_z->setVal(amrex::Real(0.0));

    amrex::MultiFab const& eps = *m_eps_mf;

    // ---- bulk anisotropic (no interface): every cell gets the SAME rows,
    //      namely the inverse of material 1. This is the production path for a
    //      homogeneous crystal, and it is what makes epsilon_tensor_* usable
    //      without an interface_function.
    if (m_has_eps_tensor && !m_has_interface) {
        auto fill = [&](amrex::MultiFab& dst, int row) {
            for (amrex::MFIter mfi(dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box const& bx = mfi.growntilebox();
                auto const& inv = dst.array(mfi);
                auto const& et1 = m_eps_tensor_mf->const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    warpx::subpixel::Mat3 e1{};
                    e1.m[0] = et1(i,j,k,0); e1.m[4] = et1(i,j,k,1); e1.m[8] = et1(i,j,k,2);
                    e1.m[1] = e1.m[3] = et1(i,j,k,3);
                    e1.m[2] = e1.m[6] = et1(i,j,k,4);
                    e1.m[5] = e1.m[7] = et1(i,j,k,5);
                    warpx::subpixel::Mat3 const chi1 = warpx::subpixel::inv3_sym(e1);
                    inv(i,j,k,0) = chi1(row,0);
                    inv(i,j,k,1) = chi1(row,1);
                    inv(i,j,k,2) = chi1(row,2);
                });
            }
        };
        fill(*m_eps_inv_x, 0);
        fill(*m_eps_inv_y, 1);
        fill(*m_eps_inv_z, 2);
        m_subpixel = true;
        if (m_has_eps2_tensor) {
            amrex::Print() << "[SUBPIXEL] bulk tensor mode: epsilon2_tensor is ignored "

                           << "(no interface_function given)\n";
        }
        return;
    }

    if (!m_has_interface) {
        // No dielectric interface: the tensor is diagonal, 1/eps.
        // Row `row` of chi^-1: the diagonal entry is 1/eps and the other two are zero.
        // (Filling only component 0 leaves chi^-1_yy and chi^-1_zz at zero, which makes
        //  Ey/Ez evolve as alpha*E and decay to nothing -- this was a real bug: it hit
        //  every uniform-medium run with subpixel=1 and no interface.)
        auto fill_diag = [&](amrex::MultiFab& dst, int row) {
            for (amrex::MFIter mfi(dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box const& bx = mfi.growntilebox();
                auto const& inv = dst.array(mfi);
                auto const& e   = eps.const_array(mfi);
                amrex::Dim3 const lo = amrex::lbound(e);
                amrex::Dim3 const hi = amrex::ubound(e);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    int const ii = amrex::min(amrex::max(i, lo.x), hi.x);
                    int const jj = amrex::min(amrex::max(j, lo.y), hi.y);
                    amrex::Real const d = amrex::Real(1.0) / e(ii, jj, k);
                    for (int c = 0; c < 3; ++c) {
                        inv(i,j,k,c) = (c == row) ? d : amrex::Real(0.0);
                    }
                });
            }
        };
        fill_diag(*m_eps_inv_x, 0);
        fill_diag(*m_eps_inv_y, 1);
        fill_diag(*m_eps_inv_z, 2);
        m_subpixel = true;
        return;
    }

    // Kottke 2008 isotropic subpixel, per cell:
    //   f      = fraction of material 1 inside the cell (sub-sampled sign of phi)
    //   n      = normalised gradient of phi (interface normal)
    //   <eps>  = f*eps1 + (1-f)*eps2
    //   <1/eps>= f/eps1 + (1-f)/eps2
    //   chi^-1 = n n^T (<1/eps> - 1/<eps>) + I/<eps>      (the 3x3 inverse tensor)
    // The three E components take the corresponding row of chi^-1 (collocated
    // approximation, same as Meep).
    // NOTE: eps1 here is only the default scalar, used for the diagnostic print and as a
    // fallback. The Kottke rows read the LOCAL epsilon from m_eps_mf(i,j,k) inside the
    // kernel, so a spatially varying epsilon (epsilon_function) is honoured. eps2 has no
    // parser in the current schema (constant only, asserted below).
    amrex::Real const eps1 = m_epsilon;
    amrex::Real const eps2 = (m_epsilon2 > 0.0) ? m_epsilon2 : m_epsilon;
    // interface_function needs a material-2 permittivity; without it the Kottke rows
    // and fill_sigma_mats would silently fall back to eps1 (wrong by O(1)).
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_has_interface
        || (m_epsilon2 > amrex::Real(0.0) || m_has_eps2_tensor),
        "macroscopic: interface_function requires a material-2 permittivity "
        "(epsilon2 or epsilon2_tensor_*)");
    auto const dx = geom.CellSize(0);
    // z is the last mesh direction: index AMREX_SPACEDIM-1, not 2 (invalid in 2D)
    auto const dz = geom.CellSize(AMREX_SPACEDIM-1);
    auto const plo = geom.ProbDomain().lo();
    auto const phi_exec = m_interface_parser->compile<3>();

    constexpr int nsub = 4;   // 4x4 sub-sampling for the fill fraction
    amrex::Real const inv_nsub = amrex::Real(1.0) / static_cast<amrex::Real>(nsub);

    bool const has_tensor = m_has_eps_tensor;
    bool const has_tensor2 = m_has_eps2_tensor;

    auto fill_rows = [&](amrex::MultiFab& dst, int row) {
        for (amrex::MFIter mfi(dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.growntilebox();
            auto const& inv = dst.array(mfi);
            auto const& ea  = m_eps_mf->const_array(mfi);
            auto const& et1 = has_tensor  ? m_eps_tensor_mf->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            auto const& et2 = has_tensor2 ? m_eps2_tensor_mf->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::Real const xc = plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx;
                // In 2D the second index is z, so the last valid component of ProbLo is z.
                amrex::Real const zc = plo[AMREX_SPACEDIM-1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dz;

                // fill fraction of material 1
                int n1 = 0;
                for (int a = 0; a < nsub; ++a) {
                    for (int b = 0; b < nsub; ++b) {
                        amrex::Real const x = xc + (static_cast<amrex::Real>(a) + amrex::Real(0.5)) * inv_nsub * dx - amrex::Real(0.5) * dx;
                        amrex::Real const z = zc + (static_cast<amrex::Real>(b) + amrex::Real(0.5)) * inv_nsub * dz - amrex::Real(0.5) * dz;
                        if (phi_exec(x, amrex::Real(0.0), z) > amrex::Real(0.0)) { ++n1; }
                    }
                }
                amrex::Real const f = static_cast<amrex::Real>(n1) * inv_nsub * inv_nsub;

                // interface normal from central differences of phi
                amrex::Real const h = amrex::Real(0.5) * amrex::min(dx, dz);
                amrex::Real const gx = phi_exec(xc + h, amrex::Real(0.0), zc) - phi_exec(xc - h, amrex::Real(0.0), zc);
                amrex::Real const gz = phi_exec(xc, amrex::Real(0.0), zc + h) - phi_exec(xc, amrex::Real(0.0), zc - h);
                amrex::Real const gn = std::sqrt(gx * gx + gz * gz);
                amrex::Real const nx = (gn > amrex::Real(1e-30)) ? gx / gn : amrex::Real(0.0);
                amrex::Real const nz = (gn > amrex::Real(1e-30)) ? gz / gn : amrex::Real(0.0);

                amrex::Real const nvec[3] = { nx, amrex::Real(0.0), nz };

                // Kottke subpixel: chi^-1 row `row`. The isotropic path uses the same
                // closed expression as before (bit-for-bit P0 regression); the
                // anisotropic path goes through the tau transform.
                amrex::Real chi[9];
                // Build both material tensors (a scalar material becomes diag(eps)) and
                // branch on POINTWISE isotropy. This is the single criterion: a scalar
                // input and an explicitly diagonal tensor therefore take the identical
                // closed-form path (bitwise), while any off-diagonal entry switches to
                // the tau transform.
                // e1/e2 read the LOCAL anisotropic tensor (6 comps: xx,yy,zz,xy,xz,yz)
                // at this cell, so a spatially varying epsilon_tensor_* is honoured.
                warpx::subpixel::Mat3 e1 = has_tensor
                    ? warpx::subpixel::Mat3{{et1(i,j,k,0), et1(i,j,k,3), et1(i,j,k,4),
                                              et1(i,j,k,3), et1(i,j,k,1), et1(i,j,k,5),
                                              et1(i,j,k,4), et1(i,j,k,5), et1(i,j,k,2)}}
                    : warpx::subpixel::diag3(ea(i,j,k));
                warpx::subpixel::Mat3 e2 = has_tensor2
                    ? warpx::subpixel::Mat3{{et2(i,j,k,0), et2(i,j,k,3), et2(i,j,k,4),
                                              et2(i,j,k,3), et2(i,j,k,1), et2(i,j,k,5),
                                              et2(i,j,k,4), et2(i,j,k,5), et2(i,j,k,2)}}
                    : warpx::subpixel::diag3(eps2);
                amrex::Real s1 = amrex::Real(0.0), s2 = amrex::Real(0.0);
                if (warpx::subpixel::is_isotropic(e1, s1) &&
                    warpx::subpixel::is_isotropic(e2, s2)) {
                    warpx::subpixel::chi_inv_isotropic(s1, s2, f, nvec, chi);
                } else {
                    warpx::subpixel::chi_inv_anisotropic(e1, e2, f, nvec, chi);
                }
                for (int c = 0; c < 3; ++c) { inv(i,j,k,c) = chi[3*row + c]; }
            });
        }
    };

    fill_rows(*m_eps_inv_x, 0);
    fill_rows(*m_eps_inv_y, 1);
    fill_rows(*m_eps_inv_z, 2);

    m_subpixel = true;

amrex::Print() << "[SUBPIXEL] eps1=" << eps1 << " eps2=" << eps2
                   << " chi^-1_yy: min=" << m_eps_inv_y->min(1) << " max=" << m_eps_inv_y->max(1)
                   << " chi^-1_x: min=" << m_eps_inv_x->min(0) << " max=" << m_eps_inv_x->max(0) << "\n";
}

