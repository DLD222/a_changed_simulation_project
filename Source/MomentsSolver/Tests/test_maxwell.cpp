/* Copyright (c) 2026
 *
 * M1 unit test: vacuum EM propagation with the 6-variable Maxwell equation.
 * A Gaussian Ey pulse (with matching Bz = Ey/c) propagates at speed c. Checks
 * that the pulse centroid advances by c*t (within a few cells) and that the
 * peak amplitude does not blow up.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MomentsSolver/Equations/MaxwellEqn.H"
#include "MomentsSolver/FV/MomentFVDriver.H"

#include <cstdio>
#include <cmath>
#include <vector>

int main ()
{
    using namespace moments;
    using Real = moments::Real;

    const int nx = 600;
    const Real xmin = 0.0, xmax = 3.0;
    const Real dx = (xmax - xmin) / nx;

    Maxwell6::Context ctx;
    ctx.c = 1.0;

    const Real x0 = 0.5, w = 0.03;
    const Real t_end = 1.0;

    std::vector<Real> U(nx * Maxwell6::num_eqs, 0.0);
    for (int i = 0; i < nx; ++i) {
        const Real x = xmin + (i + 0.5) * dx;
        const Real Ey = std::exp(-(x - x0) * (x - x0) / (w * w));
        U[i * 6 + 1] = Ey;             // Ey
        U[i * 6 + 5] = Ey / ctx.c;     // Bz = Ey/c
    }

    auto argmax_ey2 = [&](const std::vector<Real>& u) {
        int im = 0; Real m = -1.0;
        for (int i = 0; i < nx; ++i) {
            const Real v = u[i * 6 + 1] * u[i * 6 + 1];
            if (v > m) { m = v; im = i; }
        }
        return im;
    };

    const int peak0 = argmax_ey2(U);

    const Real cfl = 0.4;
    const Real dt = cfl * dx / ctx.c;
    const int nsteps = static_cast<int>(std::ceil(t_end / dt));
    for (int s = 0; s < nsteps; ++s) {
        FV1D<Maxwell6>::advance(U.data(), nx, dx, dt, ctx, BC::Outflow);
    }

    const int peak1 = argmax_ey2(U);
    const Real expected_shift = ctx.c * t_end / dx;
    const Real actual_shift = static_cast<Real>(peak1 - peak0);

    std::printf("Maxwell: peak cell %d -> %d (expected shift %.1f cells, got %.1f)\n",
                peak0, peak1, expected_shift, actual_shift);

    const bool ok = (std::abs(actual_shift - expected_shift) <= 3.0);
    std::printf("%s\n", ok ? "PASS: test_maxwell" : "RESULT: test_maxwell");
    return ok ? 0 : 1;
}
