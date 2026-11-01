/* Copyright (c) 2026
 *
 * M2 test: transverse EM wave dispersion in a two-fluid plasma.
 *
 * Electron + ion fluids (cold) in a periodic box. A transverse standing wave
 * (Ey only, Bz = 0) is seeded and oscillates at the cold-plasma EM frequency
 *     omega^2 = k^2 c^2 + omega_pe^2 + omega_pi^2 .
 * We measure the oscillation frequency of Ey and compare to this dispersion.
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
    const Real L = 2.0 * M_PI;                  // k = 1
    const Real dx = L / nx;
    const Real k = 2.0 * M_PI / L;

    const Real n0 = 1.0;
    const Real E0 = 1e-3;                       // small amplitude (linear)
    const Real p0 = 1e-4;                       // near-cold background pressure

    std::vector<Species5> sp(2);
    sp[0].q = -1.0; sp[0].m = 1.0;  sp[0].gamma = 5.0 / 3.0; // electron
    sp[1].q =  1.0; sp[1].m = 100.0; sp[1].gamma = 5.0 / 3.0; // ion (M=100)

    std::vector<Real> field(nx * 6, 0.0);
    std::vector<Real> Ue(nx * 5, 0.0), Ui(nx * 5, 0.0);
    for (int i = 0; i < nx; ++i) {
        const Real x = (i + 0.5) * dx;
        field[i * 6 + 1] = E0 * std::cos(k * x); // Ey (transverse), Bz = 0
        Ue[i * 5 + 0] = n0; Ue[i * 5 + 4] = p0 / (sp[0].gamma - 1.0); // cold, at rest
        Ui[i * 5 + 0] = n0; Ui[i * 5 + 4] = p0 / (sp[1].gamma - 1.0);
    }

    const Real wpe2 = n0 * 1.0 / 1.0;              // omega_pe^2 = n0 e^2/(m_e eps0) = 1
    const Real wpi2 = n0 * 1.0 / sp[1].m;          // omega_pi^2 = n0 e^2/(m_i eps0)
    const Real w2 = k * k * c * c + wpe2 + wpi2;
    const Real w_exact = std::sqrt(w2);
    const Real T_period = 2.0 * M_PI / w_exact;

    const Real dt = 0.4 * dx / c;
    const int nsteps = static_cast<int>(std::ceil(6.0 * T_period / dt));

    std::vector<Real*> Up{ Ue.data(), Ui.data() };
    const int probe = 0; // cos(0) = 1, Ey peak
    std::vector<Real> Ey_hist;
    Ey_hist.reserve(nsteps);

    for (int s = 0; s < nsteps; ++s) {
        strang_step(field.data(), Up, sp, nx, dx, dt, eps0, c, BC::Periodic);
        Ey_hist.push_back(field[probe * 6 + 1]);
    }

    int crossings = 0;
    for (std::size_t i = 1; i < Ey_hist.size(); ++i) {
        if (Ey_hist[i - 1] * Ey_hist[i] < 0.0) { ++crossings; }
    }
    const Real t_total = nsteps * dt;
    const Real w_meas = 2.0 * M_PI * crossings / (2.0 * t_total);

    std::printf("EM dispersion: omega_exact = %.6f (k^2 c^2=%g, wpe^2=%g, wpi^2=%g), measured = %.6f (rel err %.2e)\n",
                w_exact, k * k * c * c, wpe2, wpi2, w_meas, std::abs(w_meas - w_exact) / w_exact);

    const bool ok = (std::abs(w_meas - w_exact) / w_exact < 0.02);
    std::printf("%s\n", ok ? "PASS: test_em_dispersion" : "RESULT: test_em_dispersion");
    return ok ? 0 : 1;
}
