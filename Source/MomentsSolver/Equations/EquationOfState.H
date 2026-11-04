/* Copyright (c) 2026
 *
 * MomentSolver: pluggable equation of state (EOS) interface.
 *
 * The 5-moment Euler equation needs p(rho, e) and its derivatives to form the
 * flux and the sound speed.  Rather than hard-coding the ideal gas, the EOS is
 * a runtime table of function pointers (GPU-copyable, like the equation
 * vtable): users can plug an ideal gas, an isothermal law, a tabulated EOS, or
 * a parser-backed EOS without touching the solver.
 *
 * Conventions (all species share them):
 *   rho : mass density
 *   e   : specific internal energy (energy per mass)
 *   p   : pressure (energy per volume)
 *   T   : temperature (optional, used by the Braginskii friction)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_EQUATIONS_EQUATION_OF_STATE_H_
#define MOMENTSOLVER_EQUATIONS_EQUATION_OF_STATE_H_

#include "Equations/MomentsConfig.H"

#include <cmath>

namespace moments {

/**
 * \brief Runtime EOS: p(rho,e), (dp/drho)_e, (dp/de)_rho, e(rho,p), T(rho,e).
 *
 * All callbacks receive the opaque `ctx` pointer.  `valid()` checks that the
 * mandatory callbacks are present.  The struct is trivially copyable, so it
 * can be captured by value in GPU kernels.
 */
struct Eos
{
    Real (*pressure)(const void*, Real rho, Real e) = nullptr;
    Real (*dp_drho_e)(const void*, Real rho, Real e) = nullptr;
    Real (*dp_de_rho)(const void*, Real rho, Real e) = nullptr;
    Real (*energy_from_p)(const void*, Real rho, Real p) = nullptr; // optional
    Real (*temperature)(const void*, Real rho, Real e) = nullptr;   // optional
    const void* ctx = nullptr;

    MOMENTS_GPU_HOST_DEVICE bool valid () const noexcept
    {
        return pressure != nullptr && dp_drho_e != nullptr && dp_de_rho != nullptr;
    }
    MOMENTS_GPU_HOST_DEVICE Real p (Real rho, Real e) const noexcept
    {
        return pressure(ctx, rho, e);
    }
    MOMENTS_GPU_HOST_DEVICE Real e_from_p (Real rho, Real pv) const noexcept
    {
        return energy_from_p ? energy_from_p(ctx, rho, pv) : Real(0);
    }
    MOMENTS_GPU_HOST_DEVICE Real T (Real rho, Real e) const noexcept
    {
        return temperature ? temperature(ctx, rho, e) : Real(0);
    }
    /** Adiabatic sound speed squared: (dp/drho)_s = (dp/drho)_e + (p/rho^2)(dp/de)_rho. */
    MOMENTS_GPU_HOST_DEVICE Real cs2 (Real rho, Real e) const noexcept
    {
        const Real pv = p(rho, e);
        return dp_drho_e(ctx, rho, e)
             + (pv / (rho * rho)) * dp_de_rho(ctx, rho, e);
    }
};

/** Ideal gas parameters: p = (gamma-1) rho e, T = e/cv. */
struct IdealGasEos
{
    Real gamma = Real(5.0) / Real(3.0);
    Real cv    = Real(1.0); // e = cv T
};

MOMENTS_GPU_HOST_DEVICE inline Real ideal_gas_pressure (const void* c, Real rho, Real e) noexcept
{
    return (static_cast<const IdealGasEos*>(c)->gamma - Real(1)) * rho * e;
}
MOMENTS_GPU_HOST_DEVICE inline Real ideal_gas_dp_drho (const void* c, Real /*rho*/, Real e) noexcept
{
    return (static_cast<const IdealGasEos*>(c)->gamma - Real(1)) * e;
}
MOMENTS_GPU_HOST_DEVICE inline Real ideal_gas_dp_de (const void* c, Real rho, Real /*e*/) noexcept
{
    return (static_cast<const IdealGasEos*>(c)->gamma - Real(1)) * rho;
}
MOMENTS_GPU_HOST_DEVICE inline Real ideal_gas_e_from_p (const void* c, Real rho, Real p) noexcept
{
    return p / ((static_cast<const IdealGasEos*>(c)->gamma - Real(1)) * rho);
}
MOMENTS_GPU_HOST_DEVICE inline Real ideal_gas_temperature (const void* c, Real /*rho*/, Real e) noexcept
{
    return e / static_cast<const IdealGasEos*>(c)->cv;
}

/** Build an Eos that points at a user-owned IdealGasEos (must outlive the Eos). */
MOMENTS_GPU_HOST_DEVICE inline Eos make_ideal_gas (const IdealGasEos* params) noexcept
{
    Eos eos;
    eos.pressure      = &ideal_gas_pressure;
    eos.dp_drho_e     = &ideal_gas_dp_drho;
    eos.dp_de_rho     = &ideal_gas_dp_de;
    eos.energy_from_p = &ideal_gas_e_from_p;
    eos.temperature   = &ideal_gas_temperature;
    eos.ctx           = params;
    return eos;
}

/** Isothermal EOS: p = rho cs2 (sound speed fixed), e not thermodynamically tied to p. */
struct IsothermalEos
{
    Real cs2 = Real(1.0);
};

MOMENTS_GPU_HOST_DEVICE inline Real isothermal_pressure (const void* c, Real rho, Real /*e*/) noexcept
{
    return static_cast<const IsothermalEos*>(c)->cs2 * rho;
}
MOMENTS_GPU_HOST_DEVICE inline Real isothermal_dp_drho (const void* c, Real /*rho*/, Real /*e*/) noexcept
{
    return static_cast<const IsothermalEos*>(c)->cs2;
}
MOMENTS_GPU_HOST_DEVICE inline Real isothermal_dp_de (const void* /*c*/, Real /*rho*/, Real /*e*/) noexcept
{
    return Real(0);
}
MOMENTS_GPU_HOST_DEVICE inline Eos make_isothermal (const IsothermalEos* params) noexcept
{
    Eos eos;
    eos.pressure  = &isothermal_pressure;
    eos.dp_drho_e = &isothermal_dp_drho;
    eos.dp_de_rho = &isothermal_dp_de;
    eos.ctx       = params;
    return eos;
}

} // namespace moments

#endif // MOMENTSOLVER_EQUATIONS_EQUATION_OF_STATE_H_
