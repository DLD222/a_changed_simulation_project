/* Copyright (c) 2026
 *
 * Does the deposited current preserve the discrete Gauss law?
 *
 * This is the intermediate step between "charge-conserving deposition" and
 * "coupling to a real field solver".  The geometry is a Yee-like 1D chain
 * embedded in a 2D (nx x 2) periodic grid:
 *
 *     rho      : cell centre  (i)
 *     Ex       : x-face       (i-1/2), stored as (nx+1) values
 *     Gauss    : (Ex(i+1) - Ex(i))/dx = rho(i)/eps0
 *     Ampere   : Ex(i) -= dt * J(i)/eps0        (J on the same face)
 *
 * Two ways of building J are compared:
 *   (a) "rho*u"  : J = (q/m) * 0.5*(rho_L u_L + rho_R u_R)     (the naive one)
 *   (b) "flux"   : J = (q/m) * F_eff, the very mass flux that updated rho
 *
 * With (b), and with the Ampere update above, the discrete Gauss residual is
 * algebraically invariant:
 *     div E^{n+1} - rho^{n+1}/eps0 = div E^n - rho^n/eps0
 * so it must stay at roundoff forever. With (a) it drifts.
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

struct Problem
{
    int nx, ny;
    Real dx, eps0, c;
    Real q, m;
    Real bg;   //!< immobile background charge density (neutralises the mean)
};

void fill_state (std::vector<Real>& U, const Problem& p)
{
    const Real k = Real(2) * Real(M_PI);   // one wavelength over the box
    for (int i = 0; i < p.nx; ++i) {
        for (int j = 0; j < p.ny; ++j) {
            const Real x = (i + Real(0.5)) * p.dx;
            const Real rho = Real(1.0) + Real(0.01) * std::cos(k * x);
            const Real Eint = Real(1e-4) / (Real(5.0) / Real(3.0) - Real(1.0));
            const std::size_t c = static_cast<std::size_t>(i) * p.ny + j;
            U[c * N + 0] = rho;
            U[c * N + 1] = Real(0);
            U[c * N + 2] = Real(0);
            U[c * N + 3] = Real(0);
            U[c * N + 4] = Eint;
        }
    }
}

/** Ex on faces from Gauss, zero mean along x. */
std::vector<Real> init_face_E (const std::vector<Real>& U, const Problem& p)
{
    std::vector<Real> Ex(p.nx + 1, Real(0));
    for (int i = 1; i <= p.nx; ++i) {
        const std::size_t c = static_cast<std::size_t>(i - 1) * p.ny;
        const Real rho_q = (p.q / p.m) * U[c * N + 0] + p.bg;
        Ex[i] = Ex[i - 1] + p.dx * rho_q / p.eps0;
    }
    // with a neutralising background the chain closes: Ex[nx] == Ex[0] to roundoff
    return Ex;
}

Real gauss_residual (const std::vector<Real>& U, const std::vector<Real>& Ex, const Problem& p)
{
    Real worst = Real(0);
    for (int i = 0; i < p.nx; ++i) {
        const std::size_t c = static_cast<std::size_t>(i) * p.ny;
        const Real divE = (Ex[i + 1] - Ex[i]) / p.dx;
        const Real rho_over_eps = ((p.q / p.m) * U[c * N + 0] + p.bg) / p.eps0;
        worst = std::max(worst, std::abs(divE - rho_over_eps));
    }
    return worst;
}

/** Half Lorentz kick with the cell-centred E (average of the two faces). */
void lorentz_half (std::vector<Real>& U, const std::vector<Real>& Ex, const Problem& p, Real dt)
{
    for (int i = 0; i < p.nx; ++i) {
        const Real Ec = Real(0.5) * (Ex[i] + Ex[i + 1]);
        for (int j = 0; j < p.ny; ++j) {
            std::size_t c = static_cast<std::size_t>(i) * p.ny + j;
            const Real rho = U[c * N + 0];
            const Real qm = p.q / p.m;
            const Real ux0 = U[c * N + 1];
            U[c * N + 1] = ux0 + dt * qm * rho * Ec;
            U[c * N + 4] += dt * qm * rho * Real(0.5) * (ux0 / rho + U[c * N + 1] / rho) * Ec;
        }
    }
}

} // namespace

int main ()
{
    Problem p;
    p.nx = 64; p.ny = 2;
    p.dx = Real(1.0) / p.nx;
    p.eps0 = Real(1.0);
    p.c = Real(1.0);
    p.q = Real(-1.0);
    p.m = Real(1.0);
    p.bg = Real(1.0);   // immobile ions, +1 to cancel the mean electron charge

    const Real dt = Real(2e-3);
    const int nsteps = 400;
    Eq::Context ctx;
    BC bc[3] = {BC::Periodic, BC::Periodic, BC::Periodic};
    FVND<Eq, 2> fv(p.nx, p.ny, 1, p.dx, p.dx, Real(1.0));

    auto run = [&](bool conservative, Real& res0, Real& resMax, Real& contMax) {
        std::vector<Real> U(fv.size());
        fill_state(U, p);
        std::vector<Real> Ex = init_face_E(U, p);
        res0 = gauss_residual(U, Ex, p);
        resMax = res0;
        contMax = Real(0);

        for (int s = 0; s < nsteps; ++s) {
            lorentz_half(U, Ex, p, Real(0.5) * dt);

            std::vector<Real> Jface(p.nx + 1, Real(0));
            if (conservative) {
                YeeDeposit<2> dep;
                const Real cont = deposit_conservative<Eq, 2>(fv, U.data(), dt, ctx, bc, p.q, p.m, dep);
                contMax = std::max(contMax, cont);
                // face q of the chain: lower face of cell q (or upper face of q-1)
                for (int q = 0; q <= p.nx; ++q) {
                    const int cell = (q == p.nx) ? (p.nx - 1) : q;
                    const std::size_t c = static_cast<std::size_t>(cell) * p.ny;
                    Jface[q] = (q == p.nx) ? dep.Jhi[0][c] : dep.Jlo[0][c];
                }
            } else {
                // naive: rho*u averaged to the face
                std::vector<Real> rho(nsteps * 0 + p.nx * p.ny), ux(rho.size());
                for (std::size_t c = 0; c < rho.size(); ++c) {
                    rho[c] = U[c * N + 0];
                    ux[c]  = U[c * N + 1] / rho[c];
                }
                for (int q = 0; q <= p.nx; ++q) {
                    const int iL = (q == 0) ? p.nx - 1 : q - 1;
                    const int iR = (q == p.nx) ? 0 : q;
                    const std::size_t cL = static_cast<std::size_t>(iL) * p.ny;
                    const std::size_t cR = static_cast<std::size_t>(iR) * p.ny;
                    Jface[q] = Real(0.5) * (p.q / p.m) * (rho[cL] * ux[cL] + rho[cR] * ux[cR]);
                }
                // no transport in the naive branch: advance the fluid directly so the
                // two runs differ only in how J is built
                fv.advance(U.data(), dt, ctx, bc);
            }

            for (int q = 0; q <= p.nx; ++q) { Ex[q] -= dt * Jface[q] / p.eps0; }

            lorentz_half(U, Ex, p, Real(0.5) * dt);
            resMax = std::max(resMax, gauss_residual(U, Ex, p));
        }
    };

    Real r0c, rMc, cMc, r0n, rMn, cMn;
    run(true,  r0c, rMc, cMc);
    run(false, r0n, rMn, cMn);

    std::printf("initial Gauss residual      : %.3e\n", static_cast<double>(r0c));
    std::printf("conservative J  -> max|Gauss| = %.3e   max|continuity| = %.3e\n",
                static_cast<double>(rMc), static_cast<double>(cMc));
    std::printf("rho*u J         -> max|Gauss| = %.3e\n", static_cast<double>(rMn));

    const bool ok = (rMc <= Real(1e-12) * std::max(r0c, Real(1)))
                 && (rMn > Real(10) * std::max(rMc, Real(1e-16)));
    std::printf("%s\n", ok
        ? "PASS: conservative J preserves discrete Gauss; rho*u J drifts"
        : "FAIL: expected the conservative J to be exact and the naive one to drift");
    return ok ? 0 : 1;
}
