/* Copyright (c) 2026
 *
 * M1 unit test: 10-moment equilibrium preservation.
 * A uniform anisotropic-pressure state (Pxx != Pyy = Pzz, zero velocity) must
 * be an exact steady state of the hyperbolic 10-moment flux: the FV update
 * must leave it unchanged to machine precision over many steps.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MomentsSolver/Equations/TenMomentEqn.H"
#include "MomentsSolver/FV/MomentFVDriver.H"

#include <cstdio>
#include <cmath>
#include <vector>

int main ()
{
    using namespace moments;
    using Real = moments::Real;

    const int nx = 64;
    const Real dx = 1.0;
    const Real dt = 0.05;

    TenMoment::Context ctx;
    ctx.k0 = 1.0;

    std::vector<Real> U(nx * TenMoment::num_eqs);
    for (int i = 0; i < nx; ++i) {
        // rho=1, u=0, Pxx=2, Pyy=1, Pzz=1 (others 0)
        Real* u = &U[i * 10];
        u[0] = 1.0;            // rho
        u[1] = u[2] = u[3] = 0.0; // momentum
        u[4] = 2.0;            // Pxx_c
        u[5] = 0.0;            // Pxy_c
        u[6] = 0.0;            // Pxz_c
        u[7] = 1.0;            // Pyy_c
        u[8] = 0.0;            // Pyz_c
        u[9] = 1.0;            // Pzz_c
    }

    std::vector<Real> U0 = U;
    for (int s = 0; s < 100; ++s) {
        FV1D<TenMoment>::advance(U.data(), nx, dx, dt, ctx, BC::Periodic);
    }

    Real max_err = 0.0;
    for (int k = 0; k < nx * TenMoment::num_eqs; ++k) {
        max_err = std::max(max_err, std::abs(U[k] - U0[k]));
    }

    std::printf("TenMoment equilibrium: max change over 100 steps = %.3e\n", max_err);
    const bool ok = (max_err < 1e-12);
    std::printf("%s\n", ok ? "PASS: test_ten_moment" : "RESULT: test_ten_moment");
    return ok ? 0 : 1;
}
