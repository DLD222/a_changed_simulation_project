/* Copyright (c) 2026
 *
 * Equation-agnostic 1D finite-volume driver for the MomentSolver.
 *
 * Scheme: MUSCL-Hancock-free reconstruction (MC limiter) + local Lax-Friedrichs
 * (LLF) numerical flux + strong-stability-preserving RK3 time stepping.
 * LLF needs only max_wave_speed(), which every equation object provides, so
 * this driver runs unchanged on Euler5 / TenMoment / Maxwell6 / ...
 *
 * Boundary conditions: Periodic or Outflow (zero-gradient).
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_FV_MOMENT_FV_DRIVER_H_
#define MOMENTSOLVER_FV_MOMENT_FV_DRIVER_H_

#include "Equations/MomentsConfig.H"

#include <vector>
#include <cmath>
#include <algorithm>

namespace moments {

enum class BC : int { Periodic = 0, Outflow = 1 };

template <typename Eq>
struct FV1D
{
    using Real = moments::Real;
    using Ctx  = typename Eq::Context;
    static constexpr int N = Eq::num_eqs;

    static Real minmod (Real a, Real b)
    {
        if (a * b <= Real(0.0)) { return Real(0.0); }
        return (std::abs(a) < std::abs(b)) ? a : b;
    }

    static Real mc_limiter (Real a, Real b)
    {
        // minmod( 2a, 2b, (a+b)/2 )
        return minmod(minmod(Real(2.0) * a, Real(2.0) * b), Real(0.5) * (a + b));
    }

    /** Copy cell i (with BC handling) into `out[N]`. */
    static void get_cell (const Real* U, int nx, int i, BC bc, Real* out)
    {
        int ii = i;
        if (bc == BC::Periodic) {
            ii = ((i % nx) + nx) % nx;
        } else {
            if (i < 0)  { ii = 0; }
            if (i >= nx){ ii = nx - 1; }
        }
        for (int v = 0; v < N; ++v) { out[v] = U[ii * N + v]; }
    }

    /** RHS = -div F at every cell, into out[nx*N]. */
    static void rhs (const Real* U, int nx, Real dx, const Ctx& ctx, BC bc, Real* out)
    {
        // limited slopes
        std::vector<Real> slope(nx * N, Real(0));
        for (int i = 0; i < nx; ++i) {
            Real Um[N], U0[N], Up[N];
            get_cell(U, nx, i - 1, bc, Um);
            get_cell(U, nx, i,     bc, U0);
            get_cell(U, nx, i + 1, bc, Up);
            for (int v = 0; v < N; ++v) {
                slope[i * N + v] = mc_limiter(U0[v] - Um[v], Up[v] - U0[v]);
            }
        }

        // interface fluxes (nx+1 interfaces)
        std::vector<Real> F((nx + 1) * N, Real(0));
        for (int j = 0; j <= nx; ++j) {
            Real Uim1[N], Uj[N];
            get_cell(U, nx, j - 1, bc, Uim1);
            get_cell(U, nx, j,     bc, Uj);

            Real UL[N], UR[N];
            for (int v = 0; v < N; ++v) {
                const Real sL = (j - 1 >= 0 && j - 1 < nx) ? slope[(j - 1) * N + v] : Real(0);
                const Real sR = (j     >= 0 && j     < nx) ? slope[ j      * N + v] : Real(0);
                UL[v] = Uim1[v] + Real(0.5) * sL;
                UR[v] = Uj[v]   - Real(0.5) * sR;
            }

            Eq::numerical_flux(UL, UR, 0, F.data() + j * N, ctx);
        }

        for (int i = 0; i < nx; ++i) {
            for (int v = 0; v < N; ++v) {
                out[i * N + v] = -(F[(i + 1) * N + v] - F[i * N + v]) / dx;
            }
        }
    }

    /** Advance U[nx*N] by one dt using SSP-RK3. */
    static void advance (Real* U, int nx, Real dx, Real dt, const Ctx& ctx, BC bc)
    {
        const int sz = nx * N;
        std::vector<Real> U1(sz), U2(sz), K(sz), Un(sz);
        for (int k = 0; k < sz; ++k) { Un[k] = U[k]; }

        // stage 1
        rhs(Un.data(), nx, dx, ctx, bc, K.data());
        for (int k = 0; k < sz; ++k) { U1[k] = Un[k] + dt * K[k]; }
        // stage 2
        rhs(U1.data(), nx, dx, ctx, bc, K.data());
        for (int k = 0; k < sz; ++k) { U2[k] = Real(0.75) * Un[k] + Real(0.25) * (U1[k] + dt * K[k]); }
        // stage 3
        rhs(U2.data(), nx, dx, ctx, bc, K.data());
        for (int k = 0; k < sz; ++k) { U[k] = (Real(1) / Real(3)) * Un[k] + (Real(2) / Real(3)) * (U2[k] + dt * K[k]); }
    }
};

} // namespace moments

#endif // MOMENTSOLVER_FV_MOMENT_FV_DRIVER_H_
