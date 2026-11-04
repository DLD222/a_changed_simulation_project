/* Copyright (c) 2026
 *
 * Charge-conserving deposition tests.
 *
 *   1. advance_with_mass_flux must reproduce advance bit for bit (the RK3
 *      arithmetic is the same; only the flux is additionally recorded).
 *   2. The discrete continuity identity
 *          (rho^{n+1} - rho^n)/dt + div J = 0
 *      must hold to machine precision, because J is built from the very mass
 *      flux that updated rho.
 *   3. The deposited current must obey the 90-degree rotation symmetry.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FV/MomentFVNDDriver.H"
#include "Coupling/YeeDeposition.H"
#include "Equations/Euler5Eqn.H"

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace moments;
using Real = moments::Real;
using Eq   = Euler5;

namespace {

constexpr int N = Eq::num_eqs;

/** A smooth, fully 2D state: no symmetry that could hide an axis bug. */
void fill_state (std::vector<Real>& U, int nx, int ny, Real dx, Real dy)
{
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const Real x = (i + Real(0.5)) * dx;
            const Real y = (j + Real(0.5)) * dy;
            const Real rho = Real(1.0) + Real(0.2) * std::sin(Real(2) * Real(M_PI) * x)
                                     * std::cos(Real(2) * Real(M_PI) * y);
            const Real ux = Real(0.3) * std::cos(Real(2) * Real(M_PI) * y);
            const Real uy = Real(0.2) * std::sin(Real(2) * Real(M_PI) * x);
            const Real p  = Real(1.0) + Real(0.1) * std::sin(Real(2) * Real(M_PI) * (x + y));
            const Real E  = p / (Real(5.0) / Real(3.0) - Real(1.0))
                          + Real(0.5) * rho * (ux * ux + uy * uy);
            const std::size_t c = static_cast<std::size_t>(i) * ny + j;
            U[c * N + 0] = rho;
            U[c * N + 1] = rho * ux;
            U[c * N + 2] = rho * uy;
            U[c * N + 3] = Real(0);
            U[c * N + 4] = E;
        }
    }
}

} // namespace

int main ()
{
    const int nx = 48, ny = 48;
    const Real dx = Real(1.0) / nx, dy = Real(1.0) / ny;
    const Real dt = Real(1e-3);
    const Real q = Real(-1.0), m = Real(1.0);

    Eq::Context ctx;
    BC bc[3] = {BC::Periodic, BC::Periodic, BC::Periodic};

    // ---- 1. advance_with_mass_flux == advance ----
    {
        FVND<Eq, 2> fv(nx, ny, 1, dx, dy, Real(1.0));
        std::vector<Real> Ua(fv.size()), Ub(fv.size());
        fill_state(Ua, nx, ny, dx, dy);
        Ub = Ua;

        fv.advance(Ua.data(), dt, ctx, bc);
        YeeDeposit<2> dep;
        deposit_conservative<Eq, 2>(fv, Ub.data(), dt, ctx, bc, q, m, dep);

        Real worst = Real(0);
        for (std::size_t c = 0; c < Ua.size(); ++c) {
            worst = std::max(worst, std::abs(Ua[c] - Ub[c]));
        }
        std::printf("advance vs advance_with_mass_flux: max|dU| = %.3e\n",
                    static_cast<double>(worst));
        std::printf("%s\n", (worst == Real(0)) ? "PASS: identical state update"
                                              : "FAIL: state differs");
    }

    // ---- 2. discrete continuity ----
    {
        FVND<Eq, 2> fv(nx, ny, 1, dx, dy, Real(1.0));
        std::vector<Real> U(fv.size());
        fill_state(U, nx, ny, dx, dy);
        YeeDeposit<2> dep;
        const Real res = deposit_conservative<Eq, 2>(fv, U.data(), dt, ctx, bc, q, m, dep);

        // scale: |d rho/dt| is O(1) here, so the residual must be at roundoff
        Real ref = Real(0);
        for (std::size_t c = 0; c < dep.rho_n.size(); ++c) {
            ref = std::max(ref, std::abs(dep.rho_np1[c] - dep.rho_n[c]) / dt);
        }
        std::printf("continuity: max|(rho^{n+1}-rho^n)/dt + div J| = %.3e"
                    "   (scale |drho/dt| = %.3e)\n",
                    static_cast<double>(res), static_cast<double>(ref));
        std::printf("%s\n", (res <= Real(1e-12) * std::max(ref, Real(1)))
            ? "PASS: charge continuity exact to roundoff"
            : "FAIL: continuity violated");
    }

    // ---- 3. rotation symmetry of the deposited current ----
    {
        const int n = 32;
        const Real h = Real(1.0) / n;
        BC bcp[3] = {BC::Periodic, BC::Periodic, BC::Periodic};

        // state symmetric under x<->y (so the deposition must be too)
        auto fill_sym = [&](std::vector<Real>& U) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    const Real x = (i + Real(0.5)) * h;
                    const Real y = (j + Real(0.5)) * h;
                    const Real s = std::sin(Real(2) * Real(M_PI) * x)
                                 * std::sin(Real(2) * Real(M_PI) * y);
                    const Real rho = Real(1.0) + Real(0.2) * s;
                    const Real ux = Real(0.3) * s;
                    const Real uy = Real(0.3) * s;
                    const Real p  = Real(1.0) + Real(0.1) * s;
                    const Real E  = p / (Real(5.0) / Real(3.0) - Real(1.0))
                                  + Real(0.5) * rho * (ux * ux + uy * uy);
                    const std::size_t c = static_cast<std::size_t>(i) * n + j;
                    U[c * N + 0] = rho; U[c * N + 1] = rho * ux;
                    U[c * N + 2] = rho * uy; U[c * N + 3] = Real(0);
                    U[c * N + 4] = E;
                }
            }
        };

        FVND<Eq, 2> fv(n, n, 1, h, h, Real(1.0));
        std::vector<Real> U(fv.size());
        fill_sym(U);
        YeeDeposit<2> dep;
        deposit_conservative<Eq, 2>(fv, U.data(), dt, ctx, bcp, q, m, dep);

        Real worst = Real(0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                const std::size_t a = static_cast<std::size_t>(i) * n + j;
                const std::size_t b = static_cast<std::size_t>(j) * n + i;
                worst = std::max(worst, std::abs(dep.Jlo[0][a] - dep.Jlo[1][b]));
                worst = std::max(worst, std::abs(dep.Jlo[1][a] - dep.Jlo[0][b]));
                worst = std::max(worst, std::abs(dep.rho_np1[a] - dep.rho_np1[b]));
            }
        }
        std::printf("deposited current, x<->y symmetry: max|diff| = %.3e\n",
                    static_cast<double>(worst));
        std::printf("%s\n", (worst < Real(1e-12)) ? "PASS: deposited current symmetric"
                                                 : "FAIL: asymmetry in deposition");
    }

    // ---- 4. 3D continuity (checks the DIM=3 face bookkeeping) ----
    {
        const int m = 16;
        const Real h = Real(1.0) / m;
        FVND<Eq, 3> fv(m, m, m, h, h, h);
        BC bc3[3] = {BC::Periodic, BC::Periodic, BC::Periodic};
        std::vector<Real> U(fv.size());
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < m; ++k) {
                    const Real x = (i + Real(0.5)) * h;
                    const Real y = (j + Real(0.5)) * h;
                    const Real z = (k + Real(0.5)) * h;
                    const Real rho = Real(1.0) + Real(0.1) * std::sin(Real(2) * Real(M_PI) * x)
                        * std::cos(Real(2) * Real(M_PI) * y) * std::sin(Real(2) * Real(M_PI) * z);
                    const Real ux = Real(0.2) * std::sin(Real(2) * Real(M_PI) * y);
                    const Real uy = Real(0.2) * std::cos(Real(2) * Real(M_PI) * z);
                    const Real uz = Real(0.2) * std::sin(Real(2) * Real(M_PI) * x);
                    const Real p  = Real(1.0);
                    const Real E  = p / (Real(5.0) / Real(3.0) - Real(1.0))
                                  + Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
                    const std::size_t c =
                        (static_cast<std::size_t>(i) * m + j) * m + k;
                    U[c * N + 0] = rho; U[c * N + 1] = rho * ux;
                    U[c * N + 2] = rho * uy; U[c * N + 3] = rho * uz;
                    U[c * N + 4] = E;
                }
            }
        }
        YeeDeposit<3> dep;
        const Real res = deposit_conservative<Eq, 3>(fv, U.data(), dt, ctx, bc3, q, m, dep);
        Real ref = Real(0);
        for (std::size_t c = 0; c < dep.rho_n.size(); ++c) {
            ref = std::max(ref, std::abs(dep.rho_np1[c] - dep.rho_n[c]) / dt);
        }
        std::printf("3D continuity: max|residual| = %.3e  (scale %.3e, relative %.2e)\n",
                    static_cast<double>(res), static_cast<double>(ref),
                    static_cast<double>(res / std::max(ref, Real(1))));
        std::printf("%s\n", (res <= Real(1e-12) * std::max(ref, Real(1)))
            ? "PASS: 3D charge continuity exact to roundoff"
            : "FAIL: 3D continuity violated");
    }

    return 0;
}
