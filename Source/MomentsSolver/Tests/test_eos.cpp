/* Copyright (c) 2026
 *
 * Unit test: pluggable EOS for the 5-moment Euler equation.
 * Checks:
 *   - ideal-gas EOS reproduces the gamma law (p, cs2, Sod trajectory)
 *   - isothermal EOS gives p = cs2*rho and constant sound speed
 *   - a custom (polytrope) EOS works through the function-pointer interface
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MomentsSolver/Equations/Euler5Eqn.H"
#include "MomentsSolver/Equations/EquationOfState.H"
#include "MomentsSolver/FV/MomentFVDriver.H"

#include <cstdio>
#include <cmath>
#include <vector>

using namespace moments;

// custom polytrope: p = K rho^gamma, e = p/((gamma-1) rho)
struct Poly { Real K; Real gamma; };
static Real poly_p (const void* c, Real rho, Real e)
{
    (void)e;
    const Poly* p = static_cast<const Poly*>(c);
    return p->K * std::pow(rho, p->gamma);
}
static Real poly_dp_drho (const void* c, Real rho, Real e)
{
    (void)e;
    const Poly* p = static_cast<const Poly*>(c);
    return p->K * p->gamma * std::pow(rho, p->gamma - Real(1));
}
static Real poly_dp_de (const void*, Real, Real) { return Real(0); }

int main ()
{
    int fails = 0;
    const Real gamma = Real(1.4);

    // state: rho=1, u=0, p=1
    Real U[5] = { Real(1), Real(0), Real(0), Real(0), Real(1) / (gamma - Real(1)) };

    // --- ideal gas via gamma
    Euler5::Context cg; cg.gamma = gamma;
    const Real p_g = Euler5::pressure(U, cg);
    const Real a_g = Euler5::max_wave_speed(U, 0, cg);

    // --- ideal gas via EOS
    IdealGasEos ig; ig.gamma = gamma;
    Euler5::Context ce; ce.use_eos = true; ce.eos = make_ideal_gas(&ig);
    const Real p_e = Euler5::pressure(U, ce);
    const Real a_e = Euler5::max_wave_speed(U, 0, ce);
    std::printf("ideal: p gamma=%.12f eos=%.12f  a gamma=%.12f eos=%.12f\n",
                p_g, p_e, a_g, a_e);
    if (std::abs(p_g - p_e) > Real(1e-12)) { ++fails; }
    if (std::abs(a_g - a_e) > Real(1e-12)) { ++fails; }

    // --- isothermal: p = cs2*rho, a = cs2^0.5
    IsothermalEos iso; iso.cs2 = Real(0.5);
    Euler5::Context ci; ci.use_eos = true; ci.eos = make_isothermal(&iso);
    const Real p_i = Euler5::pressure(U, ci);
    const Real a_i = Euler5::max_wave_speed(U, 0, ci);
    std::printf("iso  : p=%.12f (ex 0.5)  a=%.12f (ex %.12f)\n",
                p_i, a_i, std::sqrt(Real(0.5)));
    if (std::abs(p_i - Real(0.5)) > Real(1e-12)) { ++fails; }
    if (std::abs(a_i - std::sqrt(Real(0.5))) > Real(1e-12)) { ++fails; }

    // --- custom polytrope with gamma=1.4, K=1 (same p as the ideal gas at rho=1)
    Poly poly{ Real(1), gamma };
    Eos eos_poly;
    eos_poly.pressure = &poly_p;
    eos_poly.dp_drho_e = &poly_dp_drho;
    eos_poly.dp_de_rho = &poly_dp_de;
    eos_poly.ctx = &poly;
    Euler5::Context cp; cp.use_eos = true; cp.eos = eos_poly;
    const Real p_p = Euler5::pressure(U, cp);
    std::printf("poly : p=%.12f (ex 1)  cs2=%.12f\n", p_p, eos_poly.cs2(Real(1), Real(2.5)));
    if (std::abs(p_p - Real(1)) > Real(1e-12)) { ++fails; }

    // --- Sod with gamma vs with the ideal-gas EOS must agree
    const int nx = 200;
    const Real dx = Real(1) / nx, t_end = Real(0.15);
    std::vector<Real> Ug(nx*5), Ue(nx*5);
    for (int i = 0; i < nx; ++i) {
        const bool left = (i + 0.5) * dx < Real(0.5);
        const Real rho = left ? Real(1) : Real(0.125);
        const Real p   = left ? Real(1) : Real(0.1);
        const Real E   = p / (gamma - Real(1));
        for (int v = 0; v < 5; ++v) { Ug[i*5+v] = 0; Ue[i*5+v] = 0; }
        Ug[i*5+0] = rho; Ug[i*5+4] = E;
        Ue[i*5+0] = rho; Ue[i*5+4] = E;
    }
    const int nsteps = 20;
    const Real dt = t_end / nsteps;
    for (int n = 0; n < nsteps; ++n) {
        FV1D<Euler5>::advance(Ug.data(), nx, dx, dt, cg, BC::Outflow);
        FV1D<Euler5>::advance(Ue.data(), nx, dx, dt, ce, BC::Outflow);
    }
    Real maxdiff = 0;
    for (std::size_t i = 0; i < Ug.size(); ++i) {
        maxdiff = std::max(maxdiff, std::abs(Ug[i] - Ue[i]));
    }
    std::printf("Sod  : max|gamma - eos| = %.3e\n", maxdiff);
    if (maxdiff > Real(1e-12)) { ++fails; }

    if (fails == 0) { std::printf("EOS PASS\n"); return 0; }
    std::printf("EOS FAIL (%d)\n", fails);
    return 1;
}
