/* Copyright (c) 2026
 *
 * MomentSolver floating-point type & GPU qualifier configuration.
 *
 * The equation objects and FV driver are intentionally AMReX-independent:
 * they build and unit-test standalone with plain C++ (Real = double). When
 * integrated into WarpX, define MOMENTS_USE_AMREX to alias Real to amrex::Real
 * and MOMENTS_GPU_HOST_DEVICE to AMREX_GPU_HOST_DEVICE (both are identical in
 * the CPU/double builds anyway).
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MOMENTSOLVER_EQUATIONS_MOMENTS_CONFIG_H_
#define MOMENTSOLVER_EQUATIONS_MOMENTS_CONFIG_H_

#if defined(MOMENTS_USE_AMREX)
#  include <AMReX_REAL.H>
#  include <AMReX_GpuQualifiers.H>
#  define MOMENTS_GPU_HOST_DEVICE AMREX_GPU_HOST_DEVICE
namespace moments { using Real = amrex::Real; }
#else
#  define MOMENTS_GPU_HOST_DEVICE
namespace moments { using Real = double; }
#endif

#endif // MOMENTSOLVER_EQUATIONS_MOMENTS_CONFIG_H_
