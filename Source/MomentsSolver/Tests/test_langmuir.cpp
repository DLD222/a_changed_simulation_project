/* Copyright (c) 2026
 *
 * M2 test: electron plasma (Langmuir) oscillation.
 *
 * A single electron fluid (q=-1, m=1) against an immobile neutralizing ion
 * background, with a small density perturbation. The electric field is
 * initialized from Gauss's law; the system should oscillate at the electron
 * plasma frequency omega_pe = sqrt(n0 e^2 / (m eps0)) = 1 (normalized units).
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Coupling/MomentCoupling.H"

#include <cstdio>
#include <cmath>
#include <vector>

int main ()
{
    using namespace moments;
    using Real = moments::Real;

    const Real eps0 = 1.0, mu0 = 1.0;
    const Real c = 1.0 / std::sqrt(eps0 * mu0); // = 1
    const int nx = 200;
    const Real L = 1.0;
    const Real dx = L / nx;

    const Real n0 = 1.0;
    const Real k = 2.0 * M_PI / L;
    const Real dn = 0.01; // density perturbation amplitude
    const Real p0 = 1e-4; // small but non-zero background pressure (near-cold)

    // electron species (immobile ion background is neutral +n0, not evolved)
    std::vector<Species5> sp(1);
    sp[0].q = -1.0; sp[0].m = 1.0; sp[0].gamma = 5.0 / 3.0;

    std::vector<Real> field(nx * 6, 0.0);
    std::vector<Real> U(nx * 5, 0.0);
    for (int i = 0; i < nx; ++i) {
        const Real x = (i + 0.5) * dx;
        const Real ne = n0 + dn * std::cos(k * x);
        U[i * 5 + 0] = ne;            // rho = ne (m=1)
        U[i * 5 + 1] = U[i * 5 + 2] = U[i * 5 + 3] = 0.0; // u = 0 (cold)
        U[i * 5 + 4] = p0 / (sp[0].gamma - 1.0);           // E_tot = p0/(gamma-1)
        // Gauss: dEx/dx = -(ne - n0)/eps0  =>  Ex = -dn sin(kx) / (eps0 k)
        field[i * 6 + 0] = -dn * std::sin(k * x) / (eps0 * k);
    }

    const Real wp = std::sqrt(n0 * 1.0 * 1.0 / (1.0 * eps0)); // omega_pe = 1
    const Real T_period = 2.0 * M_PI / wp;

    const Real dt = 0.4 * dx / c;
    const int nsteps = static_cast<int>(std::ceil(4.0 * T_period / dt));

    std::vector<Real*> Up{ U.data() };
    const int probe = nx / 4; // sin(k x) = 1 there, Ex(t) is clean
    std::vector<Real> Ex_hist;
    Ex_hist.reserve(nsteps);

    for (int s = 0; s < nsteps; ++s) {
        strang_step(field.data(), Up, sp, nx, dx, dt, eps0, c, BC::Periodic);
        Ex_hist.push_back(field[probe * 6 + 0]);
    }

    // count sign changes of Ex(t) to estimate the frequency
    int crossings = 0;
    for (std::size_t i = 1; i < Ex_hist.size(); ++i) {
        if (Ex_hist[i - 1] * Ex_hist[i] < 0.0) { ++crossings; }
    }
    const Real t_total = nsteps * dt;
    const Real f_meas = crossings / (2.0 * t_total);
    const Real f_exact = wp / (2.0 * M_PI);

    std::printf("Langmuir: omega_pe = %.6f, measured = %.6f (rel err %.2e)\n",
                wp, 2.0 * M_PI * f_meas, std::abs(f_meas - f_exact) / f_exact);

    const bool ok = (std::abs(f_meas - f_exact) / f_exact < 0.02);
    std::printf("%s\n", ok ? "PASS: test_langmuir" : "RESULT: test_langmuir");
    return ok ? 0 : 1;
}
