/* Copyright (c) 2026
 *
 * Correctness gate for FVND (2D/3D finite-volume driver).
 *
 * The 2D and 3D drivers must reproduce the 1D result **bit for bit** on a
 * problem that is uniform along the extra axes: every line then sees identical
 * interface fluxes, so the extra axes' divergence is exactly zero.
 *
 * Two checks:
 *   1. Sod shock tube to t = 0.2: max|U_2D - U_1D| and max|U_3D - U_1D| == 0
 *   2. A uniform state stays uniform (no spurious y/z flux):
 *      max|du/dt| == 0 in exact arithmetic
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FV/MomentFVDriver.H"
#include "FV/MomentFVNDDriver.H"
#include "Equations/Euler5Eqn.H"

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace moments;
using Real = moments::Real;
using Eq   = Euler5;

namespace {

/** Sod initial condition, uniform along y/z. */
void fill_sod (std::vector<Real>& U, int nx, int ny, int nz)
{
    for (int i = 0; i < nx; ++i) {
        const Real xc = (i + Real(0.5)) / nx - Real(0.5);
        const Real rho = (xc < Real(0)) ? Real(1.0)   : Real(0.125);
        const Real p   = (xc < Real(0)) ? Real(1.0)   : Real(0.1);
        const Real E   = p / (Real(5.0) / Real(3.0) - Real(1.0));
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const std::size_t c = ((static_cast<std::size_t>(i) * ny + j) * nz + k) * Eq::num_eqs;
                U[c + 0] = rho;
                U[c + 1] = Real(0);
                U[c + 2] = Real(0);
                U[c + 3] = Real(0);
                U[c + 4] = E;
            }
        }
    }
}

/** Sod profile along y (uniform in x), 2D. */
void fill_sod_y (std::vector<Real>& U, int nx, int ny)
{
    for (int j = 0; j < ny; ++j) {
        const Real xc = (j + Real(0.5)) / ny - Real(0.5);
        const Real rho = (xc < Real(0)) ? Real(1.0)   : Real(0.125);
        const Real p   = (xc < Real(0)) ? Real(1.0)   : Real(0.1);
        const Real E   = p / (Real(5.0) / Real(3.0) - Real(1.0));
        for (int i = 0; i < nx; ++i) {
            const std::size_t c = ((static_cast<std::size_t>(i) * ny + j)) * Eq::num_eqs;
            U[c + 0] = rho; U[c + 1] = Real(0); U[c + 2] = Real(0);
            U[c + 3] = Real(0); U[c + 4] = E;
        }
    }
}

/** Sod profile along z (uniform in x,y), 3D. */
void fill_sod_z (std::vector<Real>& U, int nx, int ny, int nz)
{
    for (int k = 0; k < nz; ++k) {
        const Real xc = (k + Real(0.5)) / nz - Real(0.5);
        const Real rho = (xc < Real(0)) ? Real(1.0)   : Real(0.125);
        const Real p   = (xc < Real(0)) ? Real(1.0)   : Real(0.1);
        const Real E   = p / (Real(5.0) / Real(3.0) - Real(1.0));
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                const std::size_t c =
                    ((static_cast<std::size_t>(i) * ny + j) * nz + k) * Eq::num_eqs;
                U[c + 0] = rho; U[c + 1] = Real(0); U[c + 2] = Real(0);
                U[c + 3] = Real(0); U[c + 4] = E;
            }
        }
    }
}

} // namespace

int main ()
{
    const int nx = 200;
    const Real L = Real(1.0);
    const Real dx = L / nx;
    const Real t_end = Real(0.2);

    Eq::Context ctx;                    // ideal gas, gamma = 5/3
    BC bc1[3] = {BC::Outflow, BC::Outflow, BC::Outflow};

    // ---- 1D reference ----
    std::vector<Real> U1(nx * Eq::num_eqs);
    fill_sod(U1, nx, 1, 1);
    {
        Real t = Real(0);
        int step = 0;
        while (t < t_end) {
            Real smax = Real(0);
            for (int i = 0; i < nx; ++i) {
                smax = std::max(smax, Eq::max_wave_speed(U1.data() + i * Eq::num_eqs, 0, ctx));
            }
            Real dt = Real(0.3) * dx / smax;
            if (t + dt > t_end) { dt = t_end - t; }
            FV1D<Eq>::advance(U1.data(), nx, dx, dt, ctx, BC::Outflow);
            t += dt; ++step;
        }
    }

    // ---- 2D: ny = 4, uniform along y ----
    const int ny2 = 4;
    {
        FVND<Eq, 2> fv(nx, ny2, 1, dx, Real(1.0), Real(1.0));
        BC bc[3] = {BC::Outflow, BC::Periodic, BC::Periodic};
        std::vector<Real> U2(fv.size());
        fill_sod(U2, nx, ny2, 1);

        Real t = Real(0);
        while (t < t_end) {
            Real dt = fv.max_dt(U2.data(), ctx, Real(0.3));
            if (t + dt > t_end) { dt = t_end - t; }
            fv.advance(U2.data(), dt, ctx, bc);
            t += dt;
        }
        Real worst = Real(0);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny2; ++j) {
                const std::size_t c2 = ((static_cast<std::size_t>(i) * ny2 + j) * 1) * Eq::num_eqs;
                for (int v = 0; v < Eq::num_eqs; ++v) {
                    worst = std::max(worst, std::abs(U2[c2 + v] - U1[i * Eq::num_eqs + v]));
                }
            }
        }
        std::printf("2D vs 1D (Sod, t=0.2): max|dU| = %.3e\n", static_cast<double>(worst));
        std::printf("%s\n", (worst == Real(0)) ? "PASS: 2D reduces to 1D bit-for-bit"
                                              : "FAIL: 2D differs from 1D");
    }

    // ---- 3D: ny = nz = 3, uniform along y,z ----
    const int ny3 = 3, nz3 = 3;
    {
        FVND<Eq, 3> fv(nx, ny3, nz3, dx, Real(1.0), Real(1.0));
        BC bc[3] = {BC::Outflow, BC::Periodic, BC::Periodic};
        std::vector<Real> U3(fv.size());
        fill_sod(U3, nx, ny3, nz3);

        Real t = Real(0);
        while (t < t_end) {
            Real dt = fv.max_dt(U3.data(), ctx, Real(0.3));
            if (t + dt > t_end) { dt = t_end - t; }
            fv.advance(U3.data(), dt, ctx, bc);
            t += dt;
        }
        Real worst = Real(0);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny3; ++j) {
                for (int k = 0; k < nz3; ++k) {
                    const std::size_t c3 =
                        ((static_cast<std::size_t>(i) * ny3 + j) * nz3 + k) * Eq::num_eqs;
                    for (int v = 0; v < Eq::num_eqs; ++v) {
                        worst = std::max(worst, std::abs(U3[c3 + v] - U1[i * Eq::num_eqs + v]));
                    }
                }
            }
        }
        std::printf("3D vs 1D (Sod, t=0.2): max|dU| = %.3e\n", static_cast<double>(worst));
        std::printf("%s\n", (worst == Real(0)) ? "PASS: 3D reduces to 1D bit-for-bit"
                                              : "FAIL: 3D differs from 1D");
    }

    // ---- uniform state must stay uniform (the y/z divergence is exactly 0) ----
    {
        FVND<Eq, 3> fv(16, 5, 5, Real(0.01), Real(0.01), Real(0.01));
        BC bc[3] = {BC::Periodic, BC::Periodic, BC::Periodic};
        std::vector<Real> U(fv.size());
        for (std::size_t c = 0; c < fv.ncell_total(); ++c) {
            U[c * Eq::num_eqs + 0] = Real(1.0);
            U[c * Eq::num_eqs + 1] = Real(0.2);   // uniform drift
            U[c * Eq::num_eqs + 2] = Real(0);
            U[c * Eq::num_eqs + 3] = Real(0);
            U[c * Eq::num_eqs + 4] = Real(2.5);
        }
        std::vector<Real> U0 = U;
        fv.advance(U.data(), Real(1e-6), ctx, bc);
        Real worst = Real(0);
        for (std::size_t c = 0; c < U.size(); ++c) {
            worst = std::max(worst, std::abs(U[c] - U0[c]));
        }
        std::printf("uniform drift state (3D): max|dU| after 1 step = %.3e\n",
                    static_cast<double>(worst));
        std::printf("%s\n", (worst < Real(1e-12)) ? "PASS: uniform state preserved"
                                                 : "FAIL: uniform state disturbed");
    }

    // ---- axis-1 flux: Sod along y must equal the 1D Sod along x ----
    {
        FVND<Eq, 2> fv(4, nx, 1, Real(1.0), dx, Real(1.0));
        BC bc[3] = {BC::Periodic, BC::Outflow, BC::Periodic};
        std::vector<Real> U(fv.size());
        fill_sod_y(U, 4, nx);
        Real t = Real(0);
        while (t < t_end) {
            Real dt = fv.max_dt(U.data(), ctx, Real(0.3));
            if (t + dt > t_end) { dt = t_end - t; }
            fv.advance(U.data(), dt, ctx, bc);
            t += dt;
        }
        Real worst = Real(0);
        for (int j = 0; j < nx; ++j) {
            // momentum lives in component 2 for the y-swept run: map (0,2,1,3,4) -> (0,1,2,3,4)
            static const int mapY[5] = {0, 2, 1, 3, 4};
            const std::size_t c = (static_cast<std::size_t>(0) * nx + j) * Eq::num_eqs;
            for (int v = 0; v < Eq::num_eqs; ++v) {
                worst = std::max(worst, std::abs(U[c + mapY[v]] - U1[j * Eq::num_eqs + v]));
            }
        }
        std::printf("Sod along y (2D, dir=1) vs 1D Sod along x: max|dU| = %.3e\n",
                    static_cast<double>(worst));
        std::printf("%s\n", (worst == Real(0)) ? "PASS: dir=1 flux matches the 1D solver"
                                              : "FAIL: dir=1 differs");
    }

    // ---- axis-2 flux: Sod along z must equal the 1D Sod along x ----
    {
        FVND<Eq, 3> fv(3, 3, nx, Real(1.0), Real(1.0), dx);
        BC bc[3] = {BC::Periodic, BC::Periodic, BC::Outflow};
        std::vector<Real> U(fv.size());
        fill_sod_z(U, 3, 3, nx);
        Real t = Real(0);
        while (t < t_end) {
            Real dt = fv.max_dt(U.data(), ctx, Real(0.3));
            if (t + dt > t_end) { dt = t_end - t; }
            fv.advance(U.data(), dt, ctx, bc);
            t += dt;
        }
        Real worst = Real(0);
        for (int k = 0; k < nx; ++k) {
            const std::size_t c = (static_cast<std::size_t>(0) * 3 + 0) * nx + k;
            static const int mapZ[5] = {0, 3, 2, 1, 4};
            const std::size_t cc = c * Eq::num_eqs;
            for (int v = 0; v < Eq::num_eqs; ++v) {
                worst = std::max(worst, std::abs(U[cc + mapZ[v]] - U1[k * Eq::num_eqs + v]));
            }
        }
        std::printf("Sod along z (3D, dir=2) vs 1D Sod along x: max|dU| = %.3e\n",
                    static_cast<double>(worst));
        std::printf("%s\n", (worst == Real(0)) ? "PASS: dir=2 flux matches the 1D solver"
                                              : "FAIL: dir=2 differs");
    }

    // ---- genuine 2D check: a circular blast must be transpose-symmetric ----
    {
        const int m = 64;
        const Real h = Real(1.0) / m;
        FVND<Eq, 2> fv(m, m, 1, h, h, Real(1.0));
        BC bc[3] = {BC::Outflow, BC::Outflow, BC::Periodic};
        std::vector<Real> U(fv.size());
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                const Real x = (i + Real(0.5)) * h - Real(0.5);
                const Real y = (j + Real(0.5)) * h - Real(0.5);
                const bool hot = (x*x + y*y) < Real(0.05) * Real(0.05);
                const Real rho = hot ? Real(1.0) : Real(0.125);
                const Real p   = hot ? Real(1.0) : Real(0.1);
                const Real E   = p / (Real(5.0) / Real(3.0) - Real(1.0));
                const std::size_t c = fv.cidx(i, j, 0);
                U[c * Eq::num_eqs + 0] = rho;
                U[c * Eq::num_eqs + 1] = Real(0);
                U[c * Eq::num_eqs + 2] = Real(0);
                U[c * Eq::num_eqs + 3] = Real(0);
                U[c * Eq::num_eqs + 4] = E;
            }
        }
        // a few steps (~3 CFL) are enough to break a wrongly implemented axis
        for (int nstep = 0; nstep < 5; ++nstep) {
            const Real dt = fv.max_dt(U.data(), ctx, Real(0.3));
            fv.advance(U.data(), dt, ctx, bc);
        }
        Real worst = Real(0);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                const std::size_t a = fv.cidx(i, j, 0) * Eq::num_eqs;
                const std::size_t b = fv.cidx(j, i, 0) * Eq::num_eqs;
                // transpose swaps x<->y, so momentum components 1 and 2 swap
                worst = std::max(worst, std::abs(U[a + 0] - U[b + 0]));
                worst = std::max(worst, std::abs(U[a + 1] - U[b + 2]));
                worst = std::max(worst, std::abs(U[a + 2] - U[b + 1]));
                worst = std::max(worst, std::abs(U[a + 3] - U[b + 3]));
                worst = std::max(worst, std::abs(U[a + 4] - U[b + 4]));
            }
        }
        std::printf("2D circular blast: max|U(i,j) - U(j,i)| = %.3e\n",
                    static_cast<double>(worst));
        std::printf("%s\n", (worst < Real(1e-12)) ? "PASS: x<->y symmetry exact"
                                                  : "FAIL: 2D symmetry broken");
    }

    return 0;
}
