/* Copyright 2020 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_MACROSCOPICPROPERTIES_H_
#define WARPX_MACROSCOPICPROPERTIES_H_

#include "MacroscopicProperties_fwd.H"

#include "Utils/WarpXConst.H"

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>


/**
 * \brief This class contains the macroscopic properties of the medium needed to
 * evaluate macroscopic Maxwell equation.
 */
class MacroscopicProperties
{
public:
    MacroscopicProperties (); // constructor
    /** Read user-defined macroscopic properties. Called in constructor. */
    void ReadParameters ();

    /**
     * \brief Allocate multifabs storing macroscopic multifabs
     *
     * @param[in] ba the box array associated to the multifabs E and B
     * @param[in] dm the distribution mapping
     * @param[in] ng_EB_alloc guard cells allocated for multifabs E and B
     */
    void AllocateLevelMFs (
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        const amrex::IntVect& ng_EB_alloc,
        const amrex::IntVect& Ex_nodal,
        const amrex::IntVect& Ey_nodal,
        const amrex::IntVect& Ez_nodal );

    /**
     * \brief Initialize multifabs storing macroscopic multifabs
     *
     * @param[in] geom the geometry
     * @param[in] Ex_stag staggering of the Ex field
     * @param[in] Ey_stag staggering of the Ey field
     * @param[in] Ez_stag staggering of the Ez field
     */
    void InitData (
        const amrex::Geometry& geom,
        const amrex::IntVect& Ex_stag,
        const amrex::IntVect& Ey_stag,
        const amrex::IntVect& Ez_stag);

    /** return MultiFab, sigma (conductivity) of the medium. */
    amrex::MultiFab& getsigma_mf  () {return (*m_sigma_mf);}
    /** return MultiFab, epsilon (permittivity) of the medium. */
    amrex::MultiFab& getepsilon_mf  () {return (*m_eps_mf);}
    amrex::MultiFab& get_eps_node_mf (int c) {return (*m_eps_node[c]);}
    amrex::MultiFab& get_sigma_m1_node_mf (int c) {return (*m_m1_node[c]);}
    amrex::MultiFab& get_sigma_m2_node_mf (int c) {return (*m_m2_node[c]);}
    amrex::MultiFab& get_sigma_node_mf (int c) {return (*m_sigma_node[c]);}
    /** return MultiFab, mu (permeability) of the medium. */
    amrex::MultiFab& getmu_mf  () {return (*m_mu_mf);}

    /** Inverse-permittivity rows at the E-component locations (Meep chi1inv scheme).
     *  Component `d` of `m_eps_inv_x` holds (eps^-1)_{d,x} at the Ex staggering;
     *  the three of them form one row of the 3x3 tensor. In the isotropic case
     *  they are diagonal (1/eps), so enabling this costs the same as the scalar
     *  path but lets the kernel support anisotropic materials later.
     *  Allocated only when `algo.macroscopic_subpixel = 1`. */
    amrex::MultiFab& geteps_inv_x_mf () {return (*m_eps_inv_x);}
    amrex::MultiFab& geteps_inv_y_mf () {return (*m_eps_inv_y);}
    amrex::MultiFab& geteps_inv_z_mf () {return (*m_eps_inv_z);}
    /** true when the subpixel rows were allocated (kernel can use them) */
    bool has_subpixel_rows () const {return m_subpixel;}
    /** M1 = exp(-eps^-1 sigma dt) and M2 = (1-exp(-A dt)) A^-1 eps^-1, stored as
     *  9 components each on a cell-centred MultiFab. Filled only where sigma > 0;
     *  sigma = 0 cells fall back to the plain Yee branch (bit-for-bit). */
    amrex::MultiFab& get_sigma_m1_mf () {return (*m_sigma_m1);}
    amrex::MultiFab& get_sigma_m2_mf () {return (*m_sigma_m2);}
    bool has_sigma_matrices () const {return m_sigma_mat_ready;}
    amrex::Real get_sigma_dt () const {return m_sigma_dt;}
    /** Per-instance (not process-static) seed flag for rho_medium; a second run in the
     *  same process must re-seed. */
    bool rho_medium_seeded () const {return m_rho_medium_seeded;}
    void set_rho_medium_seeded () {m_rho_medium_seeded = true;}
    int& gauss_check_step () {return m_gauss_check_step;}
    /** Fill the three rows from the scalar eps (diagonal) and, if an EB level set is
     *  available, apply the Kottke normal/tangential averaging. Called after InitEB. */
    /** Second material permittivity. Together with `macroscopic.interface_function(x,y,z)`
     *  (phi > 0 = material 1) it defines a *dielectric* interface, which is independent
     *  of the PEC embedded boundaries. See 技术路径.md §5.1. */
    amrex::Real m_epsilon2 = 0.0;
    /** whether a dielectric interface (phi) was provided */
    bool m_has_interface = false;
    /** whether the 6-component tensor inputs were provided */
    bool m_has_eps_tensor = false;
    bool m_has_eps2_tensor = false;
    /** parser for the dielectric interface level set */
    std::string m_str_interface_function;
    std::unique_ptr<amrex::Parser> m_interface_parser;
    /** true when a dielectric interface is active */
    bool has_interface () const {return m_has_interface;}

    /** Anisotropic permittivity tensors (absolute SI, same units as `macroscopic.epsilon`).
     *  Stored as 6 independent entries (xx,yy,zz,xy,xz,yz); zero means "not given".
     *  Material 1 defaults to the scalar `macroscopic.epsilon` when all six are zero. */
    amrex::GpuArray<amrex::Real,6> m_eps_tensor{};
    amrex::GpuArray<amrex::Real,6> m_eps2_tensor{};
    bool has_eps_tensor () const {return m_has_eps_tensor;}
    bool has_eps2_tensor () const {return m_has_eps2_tensor;}
    /** Cell-centred anisotropic permittivity tensors (6 comps each), filled from the
     *  per-component parsers in InitData (constant or (x,y,z) function). */
    amrex::MultiFab& get_eps_tensor_mf () {return (*m_eps_tensor_mf);}
    amrex::MultiFab& get_eps2_tensor_mf () {return (*m_eps2_tensor_mf);}
    bool has_eps_tensor_mf () const {return m_eps_tensor_mf != nullptr;}
    bool has_eps2_tensor_mf () const {return m_eps2_tensor_mf != nullptr;}

    void InitSubpixelRows ( const amrex::Geometry& geom,
                            const amrex::EBFArrayBoxFactory* eb_fact,
                            const amrex::IntVect& Ex_stag,
                            const amrex::IntVect& Ey_stag,
                            const amrex::IntVect& Ez_stag,
                            amrex::Real dt );

    /** Initializes the Multifabs storing macroscopic properties
     *  with user-defined functions(x,y,z).
     */
    void InitializeMacroMultiFabUsingParser (amrex::MultiFab *macro_mf,
                                  amrex::ParserExecutor<3> const& macro_parser,
                                  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_lev,
                                  const amrex::RealBox& prob_domain_lev);

    /** Gpu Vector with index type of the conductivity multifab */
    amrex::GpuArray<int, 3> sigma_IndexType;
    /** Gpu Vector with index type of the permittivity multifab */
    amrex::GpuArray<int, 3> epsilon_IndexType;
    /** Gpu Vector with index type of the permeability multifab */
    amrex::GpuArray<int, 3> mu_IndexType;
    /** Gpu Vector with index type of the Ex multifab */
    amrex::GpuArray<int, 3> Ex_IndexType;
    /** Gpu Vector with index type of the Ey multifab */
    amrex::GpuArray<int, 3> Ey_IndexType;
    /** Gpu Vector with index type of the Ez multifab */
    amrex::GpuArray<int, 3> Ez_IndexType;
    /** Gpu Vector with index type of coarsening ratio with default value (1,1,1) */
    amrex::GpuArray<int, 3> macro_cr_ratio;

private:

    /** Conductivity, sigma, of the medium */
    amrex::Real m_sigma = 0.0;
    /** Permittivity, epsilon, of the medium */
    amrex::Real m_epsilon = PhysConst::epsilon_0;
    /** Permeability, mu, of the medium */
    amrex::Real m_mu = PhysConst::mu0;
    /** Multifab for m_sigma */
    std::unique_ptr<amrex::MultiFab> m_sigma_mf;
    /** Multifab for m_epsilon */
    std::unique_ptr<amrex::MultiFab> m_eps_mf;
    /** Multifab for m_mu */
    std::unique_ptr<amrex::MultiFab> m_mu_mf;
    /** Cell-centred anisotropic permittivity tensors (6 comps each) */
    std::unique_ptr<amrex::MultiFab> m_eps_tensor_mf;
    std::unique_ptr<amrex::MultiFab> m_eps2_tensor_mf;
    /** Per-component parser strings and parsers (constant or (x,y,z) function) */
    std::string m_str_eps_tensor[6];
    std::string m_str_eps2_tensor[6];
    std::unique_ptr<amrex::Parser> m_eps_tensor_parser[6];
    std::unique_ptr<amrex::Parser> m_eps2_tensor_parser[6];
    /** Inverse-permittivity rows (3 comps each) at the Ex/Ey/Ez staggers */
    std::unique_ptr<amrex::MultiFab> m_eps_inv_x;
    std::unique_ptr<amrex::MultiFab> m_eps_inv_y;
    std::unique_ptr<amrex::MultiFab> m_eps_inv_z;
    /** whether the subpixel rows are allocated */
    bool m_subpixel = false;
    std::unique_ptr<amrex::MultiFab> m_sigma_m1;
    std::unique_ptr<amrex::MultiFab> m_sigma_m2;
    bool m_sigma_mat_ready = false;
    /** Node-centred material properties, one MultiFab per E component (see
     *  MEDIUM_MULTIGRID_ROOTCAUSE_AND_FIX.md §4): every shared node then has a single
     *  well-defined value and both grid copies are updated by the same rule. */
    std::unique_ptr<amrex::MultiFab> m_eps_node[3];
    std::unique_ptr<amrex::MultiFab> m_sigma_node[3];
    /** Node-centred exponential matrices, 9 comps each, one set per E component. */
    std::unique_ptr<amrex::MultiFab> m_m1_node[3];
    std::unique_ptr<amrex::MultiFab> m_m2_node[3];
    /** dt used to build m_sigma_m1/m2 (asserted constant, see 技术路径 §7.7) */
    amrex::Real m_sigma_dt = 0.0;
    bool m_rho_medium_seeded = false;
    int m_gauss_check_step = 0;

    /** Stores initialization type for conductivity : constant or parser */
    std::string m_sigma_s = "constant";
    /** Stores initialization type for permittivity : constant or parser */
    std::string m_epsilon_s = "constant";
    /** Stores initialization type for permeability : constant or parser */
    std::string m_mu_s = "constant";

    /** string for storing parser function */
    std::string m_str_sigma_function;
    std::string m_str_epsilon_function;
    std::string m_str_mu_function;
    /** Parser Wrappers */
    std::unique_ptr<amrex::Parser> m_sigma_parser;
    std::unique_ptr<amrex::Parser> m_epsilon_parser;
    std::unique_ptr<amrex::Parser> m_mu_parser;

};

/**
 * \brief
 * This struct contains only static functions to compute the coefficients for the
 * Lax-Wendroff scheme of macroscopic Maxwells equations using
 * macroscopic properties, namely, conductivity (sigma), permittivity (epsilon).
 * Permeability of the material, mu, is used as (beta/mu) for the E-update
 * defined in MacroscopicEvolveECartesian().
 */
struct LaxWendroffAlgo {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real alpha (amrex::Real const sigma,
                              amrex::Real const epsilon,
                              amrex::Real dt) {
        using namespace amrex;
        const amrex::Real fac1 = 0.5_rt * sigma * dt / epsilon;
        const amrex::Real alpha = (1._rt - fac1)/(1._rt + fac1);
        return alpha;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real beta (amrex::Real const sigma,
                             amrex::Real const epsilon,
                             amrex::Real dt) {
        using namespace amrex;
        const amrex::Real fac1 = 0.5_rt * sigma * dt / epsilon;
        const amrex::Real beta = dt / ( epsilon * (1._rt + fac1) );
        return beta;
    }

};

/**
 * \brief
 * This struct contains only static functions to compute the coefficients for the
 * BackwardEuler scheme of macroscopic Maxwells equations using
 * macroscopic properties, namely, conductivity (sigma) and permittivity (epsilon).
 * Permeability of the material, mu, is used as (beta/mu) for the E-update
 * defined in MacroscopicEvolveECartesian().
 */
struct BackwardEulerAlgo {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real alpha (amrex::Real const sigma,
                              amrex::Real const epsilon,
                              amrex::Real dt) {
        using namespace amrex;
        const amrex::Real fac1 = sigma * dt / epsilon;
        const amrex::Real alpha = (1._rt)/(1._rt + fac1);
        return alpha;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real beta (amrex::Real const sigma,
                             amrex::Real const epsilon,
                             amrex::Real dt) {
        using namespace amrex;
        const amrex::Real fac1 = sigma * dt / epsilon;
        const amrex::Real beta = dt / ( epsilon * (1._rt + fac1) );
        return beta;
    }

};

#endif // WARPX_MACROSCOPIC_PROPERTIES_H_
