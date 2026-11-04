/* Copyright (c) 2026
 *
 * Unit test: inter-species friction (three models + exponential integration).
 * Checks:
 *   - exact velocity decay exp(-nu dt) for the simple drag
 *   - total momentum conservation
 *   - total energy conservation (kinetic loss -> internal)
 *   - user-defined model reproduces the same relaxation
 *   - Braginskii rate is finite/positive and scales as T^{-3/2}
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MomentsSolver/Coupling/FluidFriction.H"

#include <cstdio>
#include <cmath>
#include <array>
#include <vector>

using namespace moments;

static Real total_momentum (std::array<std::vector<Real>,2> const& U)
{
    return (U[0][1] + U[0][2] + U[0][3]) + (U[1][1] + U[1][2] + U[1][3]);
}
static Real total_energy (std::array<std::vector<Real>,2> const& U)
{
    return U[0][4] + U[1][4];
}

static void init_state (std::array<std::vector<Real>,2>& U, Real u0, Real p, Real gamma)
{
    for (int s = 0; s < 2; ++s) {
        U[s].assign(5, Real(0));
        U[s][0] = Real(1);                 // rho
        U[s][1] = (s == 0 ? u0 : -u0);     // rho ux
        U[s][4] = Real(0.5) * U[s][0] * (s == 0 ? u0*u0 : u0*u0) + p / (gamma - Real(1));
    }
}

static Real nu_user (int, int, Real, Real, const Real*, const Real*,
                     const FrictionParams*, void*)
{
    return Real(0.3);
}

int main ()
{
    const Real gamma = Real(5.0) / Real(3.0);
    const Real p = Real(1);
    const Real u0 = Real(0.5);
    const Real dt = Real(0.1);
    const Real nu = Real(1.0);

    IdealGasEos ig; ig.gamma = gamma;
    Eos eos0 = make_ideal_gas(&ig);
    Eos eos_arr[2] = { eos0, eos0 };
    Real q[2] = { Real(-1), Real(1) };
    Real m[2] = { Real(1), Real(1) };

    int fails = 0;

    // --- 1. simple drag: exact decay + conservation
    {
        std::array<std::vector<Real>,2> U;
        init_state(U, u0, p, gamma);
        const Real P0 = total_momentum(U), E0 = total_energy(U);
        Real* ptr[2] = { U[0].data(), U[1].data() };

        FrictionParams fp; fp.model = FrictionModel::SimpleDrag; fp.nu = nu;
        apply_friction(ptr, 2, 1, dt, fp, q, m, eos_arr);

        const Real us = U[0][1] / U[0][0], ut = U[1][1] / U[1][0];
        const Real decay = std::exp(-nu * dt);
        const Real us_ex = Real(0.5) * (u0 - u0) + Real(0.5) * (u0 - (-u0)) * decay;
        const Real ut_ex = Real(0.5) * (u0 - u0) - Real(0.5) * (u0 - (-u0)) * decay;
        const Real P1 = total_momentum(U), E1 = total_energy(U);

        const Real d_us = std::abs(us - us_ex);
        const Real d_ut = std::abs(ut - ut_ex);
        const Real d_P = std::abs(P1 - P0);
        const Real d_E = std::abs(E1 - E0);
        std::printf("simple: us=%.12f ex=%.12f  ut=%.12f ex=%.12f\n", us, us_ex, ut, ut_ex);
        std::printf("        dP=%.3e dE=%.3e\n", d_P, d_E);
        if (d_us > 1e-12 || d_ut > 1e-12 || d_P > 1e-12 || d_E > 1e-12) { ++fails; }
    }

    // --- 2. user-defined model reproduces simple drag
    {
        std::array<std::vector<Real>,2> U;
        init_state(U, u0, p, gamma);
        Real* ptr[2] = { U[0].data(), U[1].data() };
        FrictionParams fp; fp.model = FrictionModel::UserDefined; fp.nu_func = &nu_user;
        apply_friction(ptr, 2, 1, dt, fp, q, m, eos_arr);
        const Real us = U[0][1] / U[0][0];
        const Real us_ex = Real(0.5) * (u0 - u0) + Real(0.5) * (u0 - (-u0)) * std::exp(-Real(0.3) * dt);
        const Real d = std::abs(us - us_ex);
        std::printf("user  : us=%.12f ex=%.12f d=%.3e\n", us, us_ex, d);
        if (d > 1e-12) { ++fails; }
    }

    // --- 3. Braginskii: unit-independent checks
    {
        const Real lnL = Real(10);
        const Real nu_T1 = friction_detail::braginskii_nu(
            Real(1), Real(1), Real(1), Real(1), Real(1), Real(1), Real(1), lnL, Real(1));
        const Real nu_T4 = friction_detail::braginskii_nu(
            Real(1), Real(1), Real(1), Real(1), Real(4), Real(1), Real(1), lnL, Real(1));
        const Real nu_T9 = friction_detail::braginskii_nu(
            Real(1), Real(1), Real(1), Real(1), Real(9), Real(1), Real(1), lnL, Real(1));
        const Real ratio4 = nu_T4 / nu_T1, ratio9 = nu_T9 / nu_T1;
        const Real ex4 = std::pow(Real(4), Real(-1.5));
        const Real ex9 = std::pow(Real(9), Real(-1.5));
        std::printf("brag  : nu(4T)/nu(T)=%.6e (ex %.6e)  nu(9T)/nu(T)=%.6e (ex %.6e)\n",
                    ratio4, ex4, ratio9, ex9);
        if (std::abs(ratio4 - ex4) > Real(1e-12)) { ++fails; }
        if (std::abs(ratio9 - ex9) > Real(1e-12)) { ++fails; }

        // relaxation with the Braginskii branch: choose nu0 so nu*dt = 0.1
        // (use physical charges so |q|/e = 1, matching the Z=1 rate above)
        Real q_brag[2] = { -friction_detail::e_charge, friction_detail::e_charge };
        std::array<std::vector<Real>,2> U;
        init_state(U, u0, p, gamma);
        Real* ptr[2] = { U[0].data(), U[1].data() };
        const Real eint = U[0][4] / U[0][0] - Real(0.5) * u0 * u0; // = p/((gamma-1)rho)
        const Real T_state = eint;                                  // cv = 1
        const Real nu1_state = friction_detail::braginskii_nu(
            Real(1), Real(1), Real(1), Real(1), T_state, Real(1), Real(1), lnL, Real(1));
        FrictionParams fp; fp.model = FrictionModel::Braginskii;
        fp.lnLambda = lnL;
        fp.nu0 = Real(0.1) / (nu1_state * dt);
        apply_friction(ptr, 2, 1, dt, fp, q_brag, m, eos_arr);
        const Real us = U[0][1] / U[0][0];
        const Real us_ex = Real(0.5) * (u0 - (-u0)) * std::exp(-Real(0.1));
        std::printf("brag  : us=%.12f ex=%.12f\n", us, us_ex);
        if (std::abs(us - us_ex) > Real(1e-10)) { ++fails; }
    }

    if (fails == 0) { std::printf("FRICTION PASS\n"); return 0; }
    std::printf("FRICTION FAIL (%d)\n", fails);
    return 1;
}
