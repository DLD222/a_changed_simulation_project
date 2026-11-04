/* Copyright (c) 2026
 *
 * Equation-agnostic N-dimensional finite-volume driver for the MomentSolver.
 *
 * Same scheme as FV1D (MC-limited reconstruction + local Lax-Friedrichs flux +
 * SSP-RK3), generalised to 2D/3D by looping over the axes and reusing the
 * per-direction hooks every equation object already provides:
 *
 *     Eq::numerical_flux(UL, UR, dir, F, ctx)      dir = 0..DIM-1
 *     Eq::max_wave_speed(U, dir, ctx)
 *
 * No equation object changes: Euler5 / TenMoment / Maxwell6 already implement
 * `dir` for all axes.
 *
 * Storage: U[((i*ny + j)*nz + k)*N + v]  (nz = 1 when DIM == 2).
 *
 * Correctness gate: a state that is uniform along y (and z) must reproduce
 * FV1D **bit for bit**, because those axes contribute identical fluxes for
 * every line and their divergence is therefore exactly 0. Tests/test_fv_nd.cpp
 * checks this against the 1D Sod run.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_FV_MOMENT_FV_ND_DRIVER_H_
#define MOMENTSOLVER_FV_MOMENT_FV_ND_DRIVER_H_

#include "Equations/MomentsConfig.H"
#include "FV/MomentFVDriver.H"   // BC, FV1D<>::mc_limiter

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>

namespace moments {

/**
 * \brief N-dimensional FV driver (DIM = 2 or 3).
 *
 * \tparam Eq  equation object (Euler5, TenMoment, Maxwell6, ...)
 * \tparam DIM active dimensions (2 or 3)
 */
template <typename Eq, int DIM>
struct FVND
{
    static_assert(DIM == 2 || DIM == 3, "FVND supports DIM = 2 or 3");
    using Real = moments::Real;
    using Ctx  = typename Eq::Context;
    static constexpr int N = Eq::num_eqs;

    int nx, ny, nz;
    Real dx, dy, dz;

    FVND (int a_nx, int a_ny, int a_nz, Real a_dx, Real a_dy, Real a_dz)
        : nx(a_nx), ny(a_ny), nz((DIM == 2) ? 1 : a_nz),
          dx(a_dx), dy(a_dy), dz((DIM == 2) ? Real(1) : a_dz) {}

    std::size_t ncell_total () const
    {
        return static_cast<std::size_t>(nx) * ny * nz;
    }

    std::size_t size () const { return ncell_total() * N; }

    int  ncell (int dir) const { return (dir == 0) ? nx : ((dir == 1) ? ny : nz); }
    Real hcell (int dir) const { return (dir == 0) ? dx : ((dir == 1) ? dy : dz); }

    /** Flat cell index. */
    std::size_t cidx (int i, int j, int k) const
    {
        return (static_cast<std::size_t>(i) * ny + j) * nz + k;
    }

    /** Flat component index. */
    std::size_t idx (int i, int j, int k, int v) const { return cidx(i,j,k) * N + v; }

    /** Decompose a flat cell index. */
    void coords (std::size_t c, int& i, int& j, int& k) const
    {
        i = static_cast<int>(c / (static_cast<std::size_t>(ny) * nz));
        j = static_cast<int>((c / static_cast<std::size_t>(nz)) % ny);
        k = static_cast<int>(c % static_cast<std::size_t>(nz));
    }

    static int wrap (int q, int n, BC bc)
    {
        if (bc == BC::Periodic) { return ((q % n) + n) % n; }
        if (q < 0)   { return 0; }
        if (q >= n)  { return n - 1; }
        return q;
    }

    /** Copy cell (i,j,k) with BC wrapping/clamping into out[N]. */
    void get_cell (const Real* U, int i, int j, int k, const BC bc[3], Real* out) const
    {
        const int ii = wrap(i, nx, bc[0]);
        const int jj = wrap(j, ny, bc[1]);
        const int kk = (DIM == 2) ? 0 : wrap(k, nz, bc[2]);
        for (int v = 0; v < N; ++v) { out[v] = U[idx(ii, jj, kk, v)]; }
    }

    /** CFL-limited dt over all active directions. */
    Real max_dt (const Real* U, const Ctx& ctx, Real cfl) const
    {
        const std::size_t nc = ncell_total();
        Real smax = Real(0);
        for (std::size_t c = 0; c < nc; ++c) {
            for (int d = 0; d < DIM; ++d) {
                smax = std::max(smax, Eq::max_wave_speed(U + c * N, d, ctx));
            }
        }
        Real hmin = dx;
        if (DIM >= 2) { hmin = std::min(hmin, dy); }
        if (DIM == 3) { hmin = std::min(hmin, dz); }
        return (smax > Real(0)) ? cfl * hmin / smax : Real(0);
    }

    /** RHS = -div F over all active axes. */
    /**
     * \brief RHS = -div F.  Optionally records the mass flux (component 0) on
     *        the lower/upper face of every cell, per axis, for charge deposition.
     *
     * `Flo[d][c]` / `Fhi[d][c]` are the fluxes at face (i-1/2) and (i+1/2) of the
     * cell in axis d.  They are exact copies of the interface fluxes used to build
     * the divergence, so (Fhi - Flo)/h reproduces `out` for component 0 exactly.
     */
    void rhs (const Real* U, const Ctx& ctx, const BC bc[3], Real* out,
              Real* Flo[3] = nullptr, Real* Fhi[3] = nullptr) const
    {
        const std::size_t nc = ncell_total();
        std::fill(out, out + nc * N, Real(0));


        for (int d = 0; d < DIM; ++d) {
            const int n = ncell(d);
            const Real h = hcell(d);
            int dims[3] = {nx, ny, nz};

            // per-axis limited slopes: slope_d[c*N+v]
            std::vector<Real> sl(nc * N, Real(0));
            for (std::size_t c = 0; c < nc; ++c) {
                int i, j, k; coords(c, i, j, k);
                int ii[3] = {i, j, k};
                Real Um[N], U0[N], Up[N];
                const int base = ii[d];          // index along the sweep axis
                ii[d] = base - 1; get_cell(U, ii[0], ii[1], ii[2], bc, Um);
                ii[d] = base;     get_cell(U, ii[0], ii[1], ii[2], bc, U0);
                ii[d] = base + 1; get_cell(U, ii[0], ii[1], ii[2], bc, Up);
                for (int v = 0; v < N; ++v) {
                    sl[c * N + v] = FV1D<Eq>::mc_limiter(U0[v] - Um[v], Up[v] - U0[v]);
                }
            }

            // iterate over every line of cells along axis d
            const int n_other0 = dims[(d + 1) % 3];
            const int n_other1 = dims[(d + 2) % 3];
            for (int a = 0; a < n_other0; ++a) {
                for (int b = 0; b < n_other1; ++b) {
                    std::vector<Real> F((n + 1) * N, Real(0));
                    for (int q = 0; q <= n; ++q) {
                        int cL[3] = {0, 0, 0}, cR[3] = {0, 0, 0};
                        cL[d] = q - 1;
                        cR[d] = q;
                        cL[(d + 1) % 3] = a; cL[(d + 2) % 3] = b;
                        cR[(d + 1) % 3] = a; cR[(d + 2) % 3] = b;

                        Real UL[N], UR[N], Um[N], U0[N], Up[N];

                        // left state: cell q-1 + 0.5*slope
                        get_cell(U, cL[0], cL[1], cL[2], bc, Um);
                        {
                            int gi = cL[0], gj = cL[1], gk = cL[2];
                            int id = cL[d];
                            const bool in = (id >= 0 && id < n);
                            int jj[3] = {gi, gj, gk};
                            jj[d] = in ? id : wrap(id, n, bc[d]);
                            const std::size_t c = cidx(jj[0], jj[1], jj[2]);
                            for (int v = 0; v < N; ++v) {
                                UL[v] = Um[v] + Real(0.5) * (in ? sl[c * N + v] : Real(0));
                            }
                        }
                        // right state: cell q - 0.5*slope
                        get_cell(U, cR[0], cR[1], cR[2], bc, U0);
                        {
                            const bool in = (q >= 0 && q < n);
                            int jj[3] = {cR[0], cR[1], cR[2]};
                            jj[d] = in ? q : wrap(q, n, bc[d]);
                            const std::size_t c = cidx(jj[0], jj[1], jj[2]);
                            for (int v = 0; v < N; ++v) {
                                UR[v] = U0[v] - Real(0.5) * (in ? sl[c * N + v] : Real(0));
                            }
                        }
                        (void)Up;

                        Eq::numerical_flux(UL, UR, d, F.data() + q * N, ctx);
                    }

                    // divergence along this line (+ optional face-flux record)
                    for (int i = 0; i < n; ++i) {
                        int cc[3] = {0, 0, 0};
                        cc[d] = i;
                        cc[(d + 1) % 3] = a; cc[(d + 2) % 3] = b;
                        const std::size_t c = cidx(cc[0], cc[1], cc[2]);
                        for (int v = 0; v < N; ++v) {
                            out[c * N + v] -= (F[(i + 1) * N + v] - F[i * N + v]) / h;
                        }
                        if (Flo != nullptr && Flo[d] != nullptr) { Flo[d][c] = F[i * N + 0]; }
                        if (Fhi != nullptr && Fhi[d] != nullptr) { Fhi[d][c] = F[(i + 1) * N + 0]; }
                    }
                }
            }
        }
    }

    /** Advance U by one dt with SSP-RK3. */
    void advance (Real* U, Real dt, const Ctx& ctx, const BC bc[3]) const
    {
        const std::size_t sz = size();
        std::vector<Real> U1(sz), U2(sz), K(sz), Un(U, U + sz);

        rhs(Un.data(), ctx, bc, K.data());
        for (std::size_t c = 0; c < sz; ++c) { U1[c] = Un[c] + dt * K[c]; }

        rhs(U1.data(), ctx, bc, K.data());
        for (std::size_t c = 0; c < sz; ++c) {
            U2[c] = Real(0.75) * Un[c] + Real(0.25) * (U1[c] + dt * K[c]);
        }

        rhs(U2.data(), ctx, bc, K.data());
        for (std::size_t c = 0; c < sz; ++c) {
            U[c] = (Real(1) / Real(3)) * Un[c] + (Real(2) / Real(3)) * (U2[c] + dt * K[c]);
        }
    }
    /**
     * \brief SSP-RK3 advance that also returns the RK-weighted mass flux on every
     *        face, `F=1/6 F1 + 1/6 F2 + 2/3 F3`.
     *
     * Because the RK3 update is U^{n+1} = U^n + dt (1/6 K1 + 1/6 K2 + 2/3 K3) and
     * K = -div F, the returned flux satisfies the discrete continuity identity
     *      (rho^{n+1} - rho^n)/dt + div F = 0
     * to machine precision.  That is what makes charge-conserving deposition
     * possible: build J from this flux instead of from rho*u.
     */
    void advance_with_mass_flux (Real* U, Real dt, const Ctx& ctx, const BC bc[3],
                                 std::vector<Real> Flo_out[3],
                                 std::vector<Real> Fhi_out[3]) const
    {
        const std::size_t sz = size();
        const std::size_t nc = ncell_total();
        for (int d = 0; d < 3; ++d) {
            Flo_out[d].assign(nc, Real(0));
            Fhi_out[d].assign(nc, Real(0));
        }

        std::vector<Real> U1(sz), U2(sz), K(sz), Un(U, U + sz);

        // stage scratch: `rhs` writes the recorded fluxes, so it must NOT alias
        // the accumulators below (otherwise stage 2 overwrites stage 1).
        std::vector<Real> flo_buf[3], fhi_buf[3];
        Real* flo[3]; Real* fhi[3];
        for (int d = 0; d < 3; ++d) {
            if (d < DIM) {
                flo_buf[d].assign(nc, Real(0));
                fhi_buf[d].assign(nc, Real(0));
                flo[d] = flo_buf[d].data();
                fhi[d] = fhi_buf[d].data();
            } else {
                flo[d] = nullptr;
                fhi[d] = nullptr;
            }
        }

        const Real w1 = Real(1) / Real(6), w2 = Real(1) / Real(6), w3 = Real(2) / Real(3);

        // stage 1
        rhs(Un.data(), ctx, bc, K.data(), flo, fhi);
        for (int d = 0; d < DIM; ++d) {
            for (std::size_t c = 0; c < nc; ++c) {
                Flo_out[d][c] = w1 * flo[d][c];
                Fhi_out[d][c] = w1 * fhi[d][c];
            }
        }
        for (std::size_t c = 0; c < sz; ++c) { U1[c] = Un[c] + dt * K[c]; }

        // stage 2
        rhs(U1.data(), ctx, bc, K.data(), flo, fhi);
        for (int d = 0; d < DIM; ++d) {
            for (std::size_t c = 0; c < nc; ++c) {
                Flo_out[d][c] += w2 * flo[d][c];
                Fhi_out[d][c] += w2 * fhi[d][c];
            }
        }
        for (std::size_t c = 0; c < sz; ++c) {
            U2[c] = Real(0.75) * Un[c] + Real(0.25) * (U1[c] + dt * K[c]);
        }

        // stage 3
        rhs(U2.data(), ctx, bc, K.data(), flo, fhi);
        for (int d = 0; d < DIM; ++d) {
            for (std::size_t c = 0; c < nc; ++c) {
                Flo_out[d][c] += w3 * flo[d][c];
                Fhi_out[d][c] += w3 * fhi[d][c];
            }
        }
        for (std::size_t c = 0; c < sz; ++c) {
            U[c] = (Real(1) / Real(3)) * Un[c] + (Real(2) / Real(3)) * (U2[c] + dt * K[c]);
        }
    }
};

} // namespace moments

#endif // MOMENTSOLVER_FV_MOMENT_FV_ND_DRIVER_H_
