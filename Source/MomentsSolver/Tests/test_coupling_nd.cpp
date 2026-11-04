/* Copyright (c) 2026
 *
 * 2D fluid <-> Maxwell Strang coupling: reduction and symmetry gates.
 *
 * The physics of the coupling source is pointwise; only the hyperbolic
 * transport is dimension dependent. Hence:
 *
 *   1. A 2D run that is uniform along y must reproduce the 1D run bit for bit.
 *   2. A 2D run whose wave vector is along y is the 90-degree rotation of the
 *      run whose wave vector is along x: after transposing the arrays (and
 *      swapping the x/y components) the two must agree bit for bit.
 *
 * Both are exact statements about the discretisation, so the tolerance is 0.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Coupling/MomentCoupling.H"

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace moments;
using Real = moments::Real;

namespace {

constexpr int N5 = 5;
constexpr int N6 = 6;

struct Setup
{
    int nx;
    Real dx;
    Real eps0, c;
    std::vector<Species5> sp;
};

Setup make_setup (int nx)
{
    Setup s;
    s.nx = nx;
    s.dx = Real(1.0) / nx;
    s.eps0 = Real(1.0);
    s.c = Real(1.0);
    s.sp.resize(1);
    s.sp[0].q = Real(-1.0);
    s.sp[0].m = Real(1.0);
    s.sp[0].gamma = Real(5.0) / Real(3.0);
    return s;
}

/** Langmuir IC: density ripple along `axis` (0 = x, 1 = y). */
void fill_langmuir (std::vector<Real>& field, std::vector<Real>& U,
                    int nx, int ny, int nz, Real dx, Real dy, Real dz,
                    Real eps0, int axis)
{
    const Real kx = Real(2.0) * Real(M_PI);
    const Real dn = Real(0.01);
    const Real p0 = Real(1e-4);
    const Real gam = Real(5.0) / Real(3.0);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const Real coord = (axis == 0) ? (i + Real(0.5)) * dx
                                               : (j + Real(0.5)) * dy;
                (void)dz; (void)k;
                const Real ne = Real(1.0) + dn * std::cos(kx * coord);
                const std::size_t c =
                    ((static_cast<std::size_t>(i) * ny + j) * nz + k);
                U[c * N5 + 0] = ne;
                U[c * N5 + 1] = Real(0);
                U[c * N5 + 2] = Real(0);
                U[c * N5 + 3] = Real(0);
                U[c * N5 + 4] = p0 / (gam - Real(1.0));
                // Gauss: dE/dx = -(ne-1)/eps0 -> E = -dn sin(kx)/(eps0 k), along the axis
                const Real E = -dn * std::sin(kx * coord) / (eps0 * kx);
                field[c * N6 + 0] = (axis == 0) ? E : Real(0);
                field[c * N6 + 1] = (axis == 1) ? E : Real(0);
                field[c * N6 + 2] = Real(0);
                field[c * N6 + 3] = Real(0);
                field[c * N6 + 4] = Real(0);
                field[c * N6 + 5] = Real(0);
            }
        }
    }
}

} // namespace

int main ()
{
    const int nx = 64;
    const int ny2 = 4;
    const int nsteps = 400;

    Setup S = make_setup(nx);
    const Real dx = S.dx, dy = S.dx;

    // ---- 1D reference ----
    std::vector<Real> f1(nx * N6, Real(0)), U1(nx * N5, Real(0));
    fill_langmuir(f1, U1, nx, 1, 1, dx, dy, Real(1), S.eps0, 0);
    {
        std::vector<Real*> Up{U1.data()};
        for (int s = 0; s < nsteps; ++s) {
            strang_step(f1.data(), Up, S.sp, nx, dx, Real(1e-2), S.eps0, S.c, BC::Periodic);
        }
    }

    // ---- 2D, uniform along y ----
    {
        std::vector<Real> f2(nx * ny2 * N6, Real(0)), U2(nx * ny2 * N5, Real(0));
        fill_langmuir(f2, U2, nx, ny2, 1, dx, dy, Real(1), S.eps0, 0);
        BC bc[3] = {BC::Periodic, BC::Periodic, BC::Periodic};
        std::vector<Real*> Up{U2.data()};
        for (int s = 0; s < nsteps; ++s) {
            strang_step_nd<2>(f2.data(), Up, S.sp, nx, ny2, 1, dx, dy, Real(1),
                              Real(1e-2), S.eps0, S.c, bc);
        }
        Real wf = Real(0), wu = Real(0);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny2; ++j) {
                const std::size_t c = static_cast<std::size_t>(i) * ny2 + j;
                for (int v = 0; v < N6; ++v) {
                    wf = std::max(wf, std::abs(f2[c * N6 + v] - f1[i * N6 + v]));
                }
                for (int v = 0; v < N5; ++v) {
                    wu = std::max(wu, std::abs(U2[c * N5 + v] - U1[i * N5 + v]));
                }
            }
        }
        std::printf("2D (uniform y) vs 1D: max|dF| = %.3e  max|dU| = %.3e\n",
                    static_cast<double>(wf), static_cast<double>(wu));
        std::printf("%s\n", (wf == Real(0) && wu == Real(0))
            ? "PASS: 2D coupling reduces to 1D bit-for-bit"
            : "FAIL: 2D coupling differs from 1D");
    }

    // ---- 2D rotation symmetry: ripple along x vs along y ----
    {
        const int m = 48;
        BC bc[3] = {BC::Periodic, BC::Periodic, BC::Periodic};
        std::vector<Real> fa(m * m * N6, Real(0)), Ua(m * m * N5, Real(0));
        std::vector<Real> fb(m * m * N6, Real(0)), Ub(m * m * N5, Real(0));
        fill_langmuir(fa, Ua, m, m, 1, dx, dy, Real(1), S.eps0, 0);   // ripple along x
        fill_langmuir(fb, Ub, m, m, 1, dx, dy, Real(1), S.eps0, 1);   // ripple along y

        std::vector<Real*> Upa{Ua.data()}, Upb{Ub.data()};
        for (int s = 0; s < 200; ++s) {
            strang_step_nd<2>(fa.data(), Upa, S.sp, m, m, 1, dx, dy, Real(1),
                              Real(1e-2), S.eps0, S.c, bc);
            strang_step_nd<2>(fb.data(), Upb, S.sp, m, m, 1, dx, dy, Real(1),
                              Real(1e-2), S.eps0, S.c, bc);
        }
        Real worst = Real(0);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                const std::size_t a = (static_cast<std::size_t>(i) * m + j);
                const std::size_t b = (static_cast<std::size_t>(j) * m + i);
                // transpose swaps x<->y
                worst = std::max(worst, std::abs(fa[a * N6 + 0] - fb[b * N6 + 1]));
                worst = std::max(worst, std::abs(fa[a * N6 + 1] - fb[b * N6 + 0]));
                worst = std::max(worst, std::abs(fa[a * N6 + 2] - fb[b * N6 + 2]));
                worst = std::max(worst, std::abs(fa[a * N6 + 3] - fb[b * N6 + 4]));
                worst = std::max(worst, std::abs(fa[a * N6 + 4] - fb[b * N6 + 3]));
                worst = std::max(worst, std::abs(fa[a * N6 + 5] - fb[b * N6 + 5]));
                worst = std::max(worst, std::abs(Ua[a * N5 + 0] - Ub[b * N5 + 0]));
                worst = std::max(worst, std::abs(Ua[a * N5 + 1] - Ub[b * N5 + 2]));
                worst = std::max(worst, std::abs(Ua[a * N5 + 2] - Ub[b * N5 + 1]));
                worst = std::max(worst, std::abs(Ua[a * N5 + 4] - Ub[b * N5 + 4]));
            }
        }
        std::printf("2D rotation symmetry (x-ripple vs y-ripple, transposed): %.3e\n",
                    static_cast<double>(worst));
        std::printf("%s\n", (worst == Real(0)) ? "PASS: x<->y symmetry exact"
                                              : "FAIL: rotation symmetry broken");
    }

    // ---- 3D, uniform along y and z (checks the DIM=3 index path) ----
    {
        const int ny3 = 2, nz3 = 2;
        std::vector<Real> f3(nx * ny3 * nz3 * N6, Real(0));
        std::vector<Real> U3(nx * ny3 * nz3 * N5, Real(0));
        fill_langmuir(f3, U3, nx, ny3, nz3, dx, dy, Real(1), S.eps0, 0);
        BC bc[3] = {BC::Periodic, BC::Periodic, BC::Periodic};
        std::vector<Real*> Up{U3.data()};
        for (int s = 0; s < nsteps; ++s) {
            strang_step_nd<3>(f3.data(), Up, S.sp, nx, ny3, nz3, dx, dy, Real(1),
                              Real(1e-2), S.eps0, S.c, bc);
        }
        Real wf = Real(0), wu = Real(0);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny3; ++j) {
                for (int k = 0; k < nz3; ++k) {
                    const std::size_t c =
                        (static_cast<std::size_t>(i) * ny3 + j) * nz3 + k;
                    for (int v = 0; v < N6; ++v) {
                        wf = std::max(wf, std::abs(f3[c * N6 + v] - f1[i * N6 + v]));
                    }
                    for (int v = 0; v < N5; ++v) {
                        wu = std::max(wu, std::abs(U3[c * N5 + v] - U1[i * N5 + v]));
                    }
                }
            }
        }
        std::printf("3D (uniform y,z) vs 1D: max|dF| = %.3e  max|dU| = %.3e\n",
                    static_cast<double>(wf), static_cast<double>(wu));
        std::printf("%s\n", (wf == Real(0) && wu == Real(0))
            ? "PASS: 3D coupling reduces to 1D bit-for-bit"
            : "FAIL: 3D coupling differs from 1D");
    }

    return 0;
}
