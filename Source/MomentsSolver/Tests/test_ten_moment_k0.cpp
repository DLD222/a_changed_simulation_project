/* Copyright (c) 2026
 *
 * M3 test: ten-moment k0 regularization (relaxation toward isotropy).
 *
 * Start from a uniform anisotropic state (rho=1, u=0, Pxx=2, Pyy=1, Pzz=1)
 * and apply only the k0 relaxation. The scalar pressure p = 4/3 and the trace
 * must be conserved, and the anisotropic part must decay as exp(-nu t) with
 * nu = v_th * k0, v_th = sqrt(p/rho).
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Coupling/TenMomentClosure.H"

#include <cstdio>
#include <cmath>
#include <vector>

int main ()
{
    using namespace moments;
    using Real = moments::Real;

    const int nx = 1; // uniform state
    const Real k0 = 1.0;
    const Real rho = 1.0;
    const Real Pxx0 = 2.0, Pyy0 = 1.0, Pzz0 = 1.0;

    std::vector<Real> U(nx * 10, 0.0);
    U[0] = rho;
    U[4] = Pxx0; U[7] = Pyy0; U[9] = Pzz0;

    const Real p_scalar = (Pxx0 + Pyy0 + Pzz0) / 3.0;      // = 4/3
    const Real v_th = std::sqrt(p_scalar / rho);           // = sqrt(4/3)
    const Real nu = v_th * k0;                             // = sqrt(4/3)

    const Real t_total = 1.0;
    const Real dt = 0.01;
    const int nsteps = static_cast<int>(t_total / dt);
    for (int s = 0; s < nsteps; ++s) {
        apply_k0_relaxation(U.data(), nx, dt, k0);
    }

    const Real Pxx = U[4], Pyy = U[7], Pzz = U[9];
    const Real exp_nu = std::exp(-nu * t_total);

    const Real Pxx_ref = (Pxx0 - p_scalar) * exp_nu + p_scalar;
    const Real Pyy_ref = (Pyy0 - p_scalar) * exp_nu + p_scalar;
    const Real Pzz_ref = (Pzz0 - p_scalar) * exp_nu + p_scalar;

    const Real trace = Pxx + Pyy + Pzz;
    const Real err_xx = std::abs(Pxx - Pxx_ref);
    const Real err_tr = std::abs(trace - 3.0 * p_scalar);

    std::printf("k0 relaxation: P=(%.6f, %.6f, %.6f), ref=(%.6f, %.6f, %.6f), trace=%.6f (should be %.6f)\n",
                Pxx, Pyy, Pzz, Pxx_ref, Pyy_ref, Pzz_ref, trace, 3.0 * p_scalar);

    const bool ok = (err_xx < 1e-6) && (err_tr < 1e-6);
    std::printf("%s\n", ok ? "PASS: test_ten_moment_k0" : "RESULT: test_ten_moment_k0");
    return ok ? 0 : 1;
}
