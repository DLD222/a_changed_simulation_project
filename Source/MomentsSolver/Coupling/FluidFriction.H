/* Copyright (c) 2026
 *
 * MomentSolver: inter-species friction (momentum exchange) for the 5-moment
 * multi-fluid system.
 *
 * The friction term in the momentum equation of species s is
 *
 *   d(rho_s u_s)/dt |_fric = sum_{t != s} R_st ,
 *   R_st = nu_st * (rho_s rho_t/(rho_s+rho_t)) * (u_t - u_s)   (R_st = -R_ts)
 *
 * so total momentum is conserved.  For frozen nu_st the pair evolution is
 * linear and is integrated EXACTLY over dt:
 *
 *   u_cm = (rho_s u_s + rho_t u_t)/(rho_s+rho_t)
 *   du   = u_s - u_t  ->  du' = du exp(-nu_st dt)
 *   u_s' = u_cm + (rho_t/(rho_s+rho_t)) du'
 *   u_t' = u_cm - (rho_s/(rho_s+rho_t)) du'
 *
 * The lost kinetic energy 0.5*mu*|du|^2*(1-exp(-2 nu dt)) (mu = reduced mass
 * density) is transferred to the internal energy, so total energy is conserved.
 *
 * Three models are provided:
 *   - SimpleDrag : constant nu (cheapest, for validation)
 *   - Braginskii : Coulomb/Braginskii electron-ion-like nu ~ n Z^2 lnL / T^{3/2}
 *   - UserDefined: user callback (external interface, e.g. table/parser)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_COUPLING_FLUID_FRICTION_H_
#define MOMENTSOLVER_COUPLING_FLUID_FRICTION_H_

#include "Equations/MomentsConfig.H"
#include "Equations/EquationOfState.H"

#include <cmath>
#include <algorithm>

namespace moments {

enum class FrictionModel : int { None = 0, SimpleDrag = 1, Braginskii = 2, UserDefined = 3 };

/**
 * \brief Friction model parameters (trivially copyable, GPU-friendly).
 *
 * `nu_func` signature (UserDefined):
 *   nu(s, t, rho_s, rho_t, u_s[3], u_t[3], params, ctx)
 * where u_* are the velocities (not momentum densities).
 */
struct FrictionParams
{
    FrictionModel model = FrictionModel::None;
    Real nu = Real(0);        //!< SimpleDrag coefficient [1/time]
    Real lnLambda = Real(10); //!< Coulomb logarithm (Braginskii)
    Real Z = Real(1);         //!< charge state (Braginskii)
    Real nu0 = Real(1);       //!< overall scale (Braginskii; unit conversion/normalization)
    Real (*nu_func)(int s, int t, Real rho_s, Real rho_t,
                    const Real* u_s, const Real* u_t,
                    const FrictionParams* fp, void* ctx) = nullptr;
    void* ctx = nullptr;
};

namespace friction_detail {

// Physical constants (SI) for the Braginskii/Coulomb coefficient.
constexpr Real e_charge = Real(1.602176634e-19);
constexpr Real eps0     = Real(8.8541878128e-12);
constexpr Real kB       = Real(1.380649e-23);
constexpr Real pi       = Real(3.14159265358979323846);

/**
 * Coulomb/Braginskii-like momentum-exchange rate between s and t:
 *   nu_st = nu0 * (4 sqrt(2 pi)/3) n_t Z_s^2 Z_t^2 e^4 lnL
 *           / (sqrt(m_s) (4 pi eps0)^2 (k_B T_s)^{3/2})
 * (the s-t reduced form; T_s from the EOS, T = e/cv if no EOS temperature).
 */
MOMENTS_GPU_HOST_DEVICE inline Real braginskii_nu (
    Real rho_s, Real rho_t, Real m_s, Real m_t, Real T_s, Real Zs, Real Zt,
    Real lnLambda, Real nu0) noexcept
{
    if (rho_s <= Real(0) || m_s <= Real(0) || T_s <= Real(0)) { return Real(0); }
    const Real n_t = rho_t / m_t;
    const Real pref = Real(4) * std::sqrt(Real(2) * pi) / Real(3);
    const Real num  = pref * n_t * Zs * Zs * Zt * Zt
                    * (e_charge*e_charge*e_charge*e_charge) * lnLambda;
    const Real den  = std::sqrt(m_s)
                    * (Real(4) * pi * eps0) * (Real(4) * pi * eps0)
                    * std::pow(kB * T_s, Real(1.5));
    return nu0 * num / den;
}

} // namespace friction_detail

/**
 * \brief Apply pairwise friction over dt (exact exponential for frozen nu).
 *
 * \param U   per-species conserved arrays, each [rho, rho ux, rho uy, rho uz, E], nx cells
 * \param ns  number of species
 * \param nx  number of cells
 * \param dt  time step
 * \param fp  friction parameters
 * \param q   species charges (for the Braginskii Z; can be nullptr)
 * \param m   species masses
 * \param eos per-species EOS (same for all here; used for T in Braginskii)
 */
inline void apply_friction (Real** U, int ns, int nx, Real dt,
                            const FrictionParams& fp, const Real* q,
                            const Real* m, const Eos* eos)
{
    if (fp.model == FrictionModel::None || dt <= Real(0)) { return; }
    const int N = 5;
    for (int s = 0; s < ns; ++s) {
        for (int t = s + 1; t < ns; ++t) {
            for (int i = 0; i < nx; ++i) {
                Real* Us = U[s] + N * i;
                Real* Ut = U[t] + N * i;
                const Real rho_s = Us[0];
                const Real rho_t = Ut[0];
                if (rho_s <= Real(0) || rho_t <= Real(0)) { continue; }

                Real us[3] = { Us[1]/rho_s, Us[2]/rho_s, Us[3]/rho_s };
                Real ut[3] = { Ut[1]/rho_t, Ut[2]/rho_t, Ut[3]/rho_t };

                Real nu = Real(0);
                if (fp.model == FrictionModel::SimpleDrag) {
                    nu = fp.nu;
                } else if (fp.model == FrictionModel::UserDefined) {
                    if (fp.nu_func) { nu = fp.nu_func(s, t, rho_s, rho_t, us, ut, &fp, fp.ctx); }
                } else { // Braginskii
                    const Real eint = Us[4]/rho_s - Real(0.5)*(us[0]*us[0]+us[1]*us[1]+us[2]*us[2]);
                    Real T_s = eos ? eos[s].T(rho_s, eint) : Real(0);
                    if (T_s <= Real(0)) { T_s = eint; }
                    const Real Zs_eff = q ? std::abs(q[s])/friction_detail::e_charge : Real(1);
                    const Real Zt_eff = q ? std::abs(q[t])/friction_detail::e_charge : Real(1);
                    nu = friction_detail::braginskii_nu(rho_s, rho_t,
                             m ? m[s] : Real(1), m ? m[t] : Real(1), T_s,
                             Zs_eff, Zt_eff, fp.lnLambda, fp.nu0);
                }
                const Real decay = std::exp(-nu * dt);

                // old internal energies (per volume), before the momentum update
                const Real eint_s_old = Us[4] - Real(0.5) * rho_s
                                      * (us[0]*us[0] + us[1]*us[1] + us[2]*us[2]);
                const Real eint_t_old = Ut[4] - Real(0.5) * rho_t
                                      * (ut[0]*ut[0] + ut[1]*ut[1] + ut[2]*ut[2]);

                // old relative velocity for the dissipated kinetic energy
                Real du0sq = Real(0);
                for (int d = 0; d < 3; ++d) {
                    const Real d0 = us[d] - ut[d];
                    du0sq += d0*d0;
                }

                // common + decaying relative velocity (exact)
                const Real rho_tot = rho_s + rho_t;
                Real ucm[3], du[3];
                for (int d = 0; d < 3; ++d) {
                    ucm[d] = (rho_s*us[d] + rho_t*ut[d]) / rho_tot;
                    du[d]  = (us[d] - ut[d]) * decay;
                }
                for (int d = 0; d < 3; ++d) {
                    us[d] = ucm[d] + (rho_t/rho_tot) * du[d];
                    ut[d] = ucm[d] - (rho_s/rho_tot) * du[d];
                }

                // dissipated kinetic energy -> internal (pair sum conserved)
                const Real mu = rho_s * rho_t / rho_tot;
                const Real diss = Real(0.5) * mu * du0sq * (Real(1) - decay*decay);
                const Real share_s = rho_s / rho_tot;
                const Real share_t = rho_t / rho_tot;

                // momenta
                Us[1] = rho_s * us[0]; Us[2] = rho_s * us[1]; Us[3] = rho_s * us[2];
                Ut[1] = rho_t * ut[0]; Ut[2] = rho_t * ut[1]; Ut[3] = rho_t * ut[2];

                // energies: new kinetic + old internal + share of dissipation
                Us[4] = Real(0.5) * rho_s * (us[0]*us[0] + us[1]*us[1] + us[2]*us[2])
                      + eint_s_old + diss * share_s;
                Ut[4] = Real(0.5) * rho_t * (ut[0]*ut[0] + ut[1]*ut[1] + ut[2]*ut[2])
                      + eint_t_old + diss * share_t;
            }
        }
    }
}

} // namespace moments

#endif // MOMENTSOLVER_COUPLING_FLUID_FRICTION_H_
