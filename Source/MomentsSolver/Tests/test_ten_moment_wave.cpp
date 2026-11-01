/* Copyright (c) 2026
 *
 * M3 test: ten-moment longitudinal sound wave (anisotropic pressure wave).
 *
 * The 10-moment linearized sound speed is c_s = sqrt(3 p0 / rho0). We seed a
 * standing longitudinal wave (rho and Pxx perturbed in phase, u = 0) and check
 * that it oscillates at omega = c_s * k.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Equations/TenMomentEqn.H"
#include "FV/MomentFVDriver.H"

#include <cstdio>
#include <cmath>
#include <vector>

int main ()
{
    using namespace moments;
    using Real = moments::Real;

    const int nx = 200;
    const Real L = 1.0;
    const Real dx = L / nx;
    const Real k = 2.0 * M_PI / L;
    const Real rho0 = 1.0, p0 = 1.0;
    const Real A = 0.01; // density perturbation amplitude

    TenMoment::Context ctx;
    ctx.k0 = 0.0; // no closure in this test (pure hyperbolic wave)

    std::vector<Real> U(nx * 10, 0.0);
    for (int i = 0; i < nx; ++i) {
        const Real x = (i + 0.5) * dx;
        const Real rho = rho0 + A * std::cos(k * x);
        U[i * 10 + 0] = rho;
        U[i * 10 + 4] = p0 + 3.0 * p0 * (rho - rho0) / rho0; // Pxx_c (u=0)
        U[i * 10 + 7] = p0;                                   // Pyy_c
        U[i * 10 + 9] = p0;                                   // Pzz_c
    }

    const Real cs = std::sqrt(3.0 * p0 / rho0);
    const Real w_exact = cs * k;
    const Real T_period = 2.0 * M_PI / w_exact;

    const Real amax = cs;
    const Real dt = 0.4 * dx / amax;
    const int nsteps = static_cast<int>(std::ceil(4.0 * T_period / dt));

    const int probe = 0; // cos(0) = 1
    std::vector<Real> rho_hist;
    rho_hist.reserve(nsteps);
    for (int s = 0; s < nsteps; ++s) {
        FV1D<TenMoment>::advance(U.data(), nx, dx, dt, ctx, BC::Periodic);
        rho_hist.push_back(U[probe * 10 + 0] - rho0);
    }

    int crossings = 0;
    for (std::size_t i = 1; i < rho_hist.size(); ++i) {
        if (rho_hist[i - 1] * rho_hist[i] < 0.0) { ++crossings; }
    }
    const Real t_total = nsteps * dt;
    const Real w_meas = 2.0 * M_PI * crossings / (2.0 * t_total);

    std::printf("10-moment sound wave: omega_exact = %.6f (cs=%.4f), measured = %.6f (rel err %.2e)\n",
                w_exact, cs, w_meas, std::abs(w_meas - w_exact) / w_exact);

    const bool ok = (std::abs(w_meas - w_exact) / w_exact < 0.02);
    std::printf("%s\n", ok ? "PASS: test_ten_moment_wave" : "RESULT: test_ten_moment_wave");
    return ok ? 0 : 1;
}
