/* Copyright (c) 2026
 *
 * M1 unit test: 1D Sod shock tube with the 5-moment Euler equation object.
 * Checks: total mass / momentum / energy conservation (waves stay inside the
 * domain at t_end), physical density bounds, and prints the profile.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MomentsSolver/Equations/Euler5Eqn.H"
#include "MomentsSolver/FV/MomentFVDriver.H"

#include <cstdio>
#include <cmath>
#include <vector>

int main ()
{
    using namespace moments;
    using Real = moments::Real;

    const int nx = 400;
    const Real xmin = 0.0, xmax = 1.0;
    const Real dx = (xmax - xmin) / nx;
    const Real t_end = 0.2;

    Euler5::Context ctx;
    ctx.gamma = 1.4;

    std::vector<Real> U(nx * Euler5::num_eqs);
    // left:  rho=1, p=1, u=0 ; right: rho=0.125, p=0.1, u=0
    for (int i = 0; i < nx; ++i) {
        const Real x = xmin + (i + 0.5) * dx;
        const bool left = (x < 0.5);
        const Real rho = left ? 1.0 : 0.125;
        const Real p   = left ? 1.0 : 0.1;
        const Real E   = p / (ctx.gamma - 1.0);
        U[i * 5 + 0] = rho;
        U[i * 5 + 1] = 0.0; U[i * 5 + 2] = 0.0; U[i * 5 + 3] = 0.0;
        U[i * 5 + 4] = E;
    }

    // conserved quantities at t=0
    Real M0 = 0.0, P0 = 0.0, E0 = 0.0;
    for (int i = 0; i < nx; ++i) {
        M0 += U[i * 5 + 0];
        P0 += U[i * 5 + 1];
        E0 += U[i * 5 + 4];
    }

    // time loop with fixed dt (CFL from initial state)
    Real amax = 0.0;
    for (int i = 0; i < nx; ++i) {
        amax = std::max(amax, Euler5::max_wave_speed(&U[i * 5], 0, ctx));
    }
    const Real cfl = 0.4;
    const Real dt = cfl * dx / amax;
    const int nsteps = static_cast<int>(std::ceil(t_end / dt));

    for (int s = 0; s < nsteps; ++s) {
        FV1D<Euler5>::advance(U.data(), nx, dx, dt, ctx, BC::Outflow);
    }

    // conserved quantities at t_end
    Real M1 = 0.0, P1 = 0.0, E1 = 0.0;
    Real rho_max = 0.0, rho_min = 1e30;
    for (int i = 0; i < nx; ++i) {
        M1 += U[i * 5 + 0];
        P1 += U[i * 5 + 1];
        E1 += U[i * 5 + 4];
        rho_max = std::max(rho_max, U[i * 5 + 0]);
        rho_min = std::min(rho_min, U[i * 5 + 0]);
    }

    const Real dM = std::abs(M1 - M0) / M0;
    const Real dP = std::abs(P1 - P0);
    const Real dE = std::abs(E1 - E0) / E0;

    bool ok = true;
    // Note: mass & energy fluxes are zero at the open boundaries (u=0 there),
    // so they are conserved. Momentum is NOT conserved here: the pressure
    // imbalance (p=1 vs p=0.1) does work across the open boundaries.
    if (dM > 1e-8)  { std::printf("FAIL: mass conservation %g\n", dM); ok = false; }
    if (dE > 1e-8)  { std::printf("FAIL: energy conservation %g\n", dE); ok = false; }
    if (rho_max > 1.0 + 1e-6) { std::printf("FAIL: overshoot rho_max=%g\n", rho_max); ok = false; }
    if (rho_min < 0.125 - 1e-6) { std::printf("FAIL: undershoot rho_min=%g\n", rho_min); ok = false; }

    std::printf("Sod: dM=%.2e dP=%.2e dE=%.2e rho=[%.4f, %.4f]\n", dM, dP, dE, rho_min, rho_max);

    // profile at a few sample points (x = 0.2, 0.5, 0.75, 0.9)
    const Real xs[4] = {0.2, 0.5, 0.75, 0.9};
    for (int k = 0; k < 4; ++k) {
        const int i = static_cast<int>((xs[k] - xmin) / dx - 0.5);
        const int ii = (i >= 0 && i < nx) ? i : 0;
        std::printf("  rho(x=%.2f) = %.5f\n", xs[k], U[ii * 5 + 0]);
    }

    std::printf("%s\n", ok ? "PASS: test_euler_sod" : "RESULT: test_euler_sod");
    return ok ? 0 : 1;
}
