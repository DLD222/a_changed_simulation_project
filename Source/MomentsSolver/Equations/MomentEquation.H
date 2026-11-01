/* Copyright (c) 2026
 *
 * MomentSolver: multi-fluid moment (5-moment / 10-moment) + full-Maxwell
 * solver forked from WarpX. Fluid-only, C++.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_EQUATIONS_MOMENT_EQUATION_H_
#define MOMENTSOLVER_EQUATIONS_MOMENT_EQUATION_H_

#include "Equations/MomentsConfig.H"

/**
 * \namespace moments
 *
 * Equation-object abstraction for hyperbolic moment systems.
 *
 * Every hyperbolic model (5-moment Euler, 10-moment pressure-tensor, vacuum
 * Maxwell, ...) exposes the same minimal contract below. The finite-volume
 * (and, later, DG) drivers only talk to this contract, so the discretization
 * is equation-agnostic.
 *
 * Design notes (see design-plan-moment-solver.md):
 *  - Storage agnostic: implementations operate on plain Real arrays (no
 *    MultiFab dependency), so the physics is trivially GPU-able and reusable
 *    by a future modal-DG discretization.
 *  - The template FV driver (FV/MomentFVDriver.H) uses static dispatch
 *    through `Eq::Context`, `Eq::flux`, `Eq::max_wave_speed`, ... A runtime
 *    "vtable" (below) is provided for the WarpX integration layer.
 */
namespace moments {

/** Direction indices used by flux()/max_wave_speed(). */
enum class Dir : int { X = 0, Y = 1, Z = 2 };

/**
 * \brief Table of function pointers implementing one hyperbolic system.
 * Used by the WarpX runtime integration (host/device copyable), Gkeyll-style.
 */
struct EquationVTable
{
    int num_eqs = 0;    //!< conserved variables per species
    int num_waves = 0;  //!< number of waves in the Riemann solution

    MOMENTS_GPU_HOST_DEVICE
    Real (*max_wave_speed)(const void* ctx, const Real* U, int dir) noexcept = nullptr;

    MOMENTS_GPU_HOST_DEVICE
    void (*flux)(const void* ctx, const Real* U, int dir, Real* F) noexcept = nullptr;

    MOMENTS_GPU_HOST_DEVICE
    void (*conserved_to_primitive)(const void* ctx, const Real* U, Real* W) noexcept = nullptr;
};

} // namespace moments

#endif // MOMENTSOLVER_EQUATIONS_MOMENT_EQUATION_H_
