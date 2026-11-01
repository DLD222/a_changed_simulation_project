/* Copyright 2021-2025 Lorenzo Giacomel, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Enabled.H"

#ifdef AMREX_USE_EB

#include "EmbeddedBoundaryInit.H"

#include "EmbeddedBoundary/WarpXFaceInfoBox.H"
#include "Fields.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Box.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiCutFab.H>

namespace web = warpx::embedded_boundary;

namespace {

/** \brief Fraction of the unit segment where phi <= 0, given phi at both ends.
 *
 *  The zero crossing is located by linear interpolation, which is exact when phi
 *  is linear along the segment (which it is inside an EB2 cell, where the level
 *  set is linearised). See C_SIDE_LS_ECT_RESEARCH.md.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real open_fraction_1d (amrex::Real p0, amrex::Real p1)
{
    bool const o0 = (p0 <= amrex::Real(0.0));
    bool const o1 = (p1 <= amrex::Real(0.0));
    if (o0 && o1) { return amrex::Real(1.0); }
    if (!o0 && !o1) { return amrex::Real(0.0); }
    amrex::Real const t = p0 / (p0 - p1);           // crossing, in (0,1)
    return o0 ? t : (amrex::Real(1.0) - t);
}

/** \brief Clamp an index to the level-set array bounds (the dual arrays' ghost
 *  cells can extend beyond the EB factory's level-set ghosts; without this the
 *  scheme-4 geometry reads out of bounds at large grids, see the 640^2 Mie run). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clamp_idx (int v, int lo, int hi)
{
    return amrex::max(lo, amrex::min(v, hi));
}

/** \brief Open length along x on the dual line through (x_i, z_{j+1/2}).
 *
 *  phi is nodal and positive inside the solid (WarpX convention), so the open
 *  part is phi <= 0. The line spans [x_i - h/2, x_i + h/2] and crosses cells
 *  (i-1,j) and (i,j); z is averaged between nodes j and j+1 because the dual
 *  line sits at the cell centre in z.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real open_len_x_from_ls (amrex::Array4<amrex::Real const> const& phi,
                                amrex::Dim3 const phi_lo, amrex::Dim3 const phi_hi,
                                int i, int j, int k, amrex::Real h, int nseg)
{
    // The write box may carry ghost cells in the (degenerate) z direction; the
    // level set only exists at k = 0 in 2D, so clamp before indexing.
    int const kk = clamp_idx(k, phi_lo.z, phi_hi.z);
    amrex::Real len = amrex::Real(0.0);
    amrex::Real prev = amrex::Real(0.0);
    bool all_open = true;
    for (int s = 0; s <= nseg; ++s) {
        amrex::Real const u = -amrex::Real(0.5)
            + static_cast<amrex::Real>(s) / static_cast<amrex::Real>(nseg);
        int   const im0 = (u < amrex::Real(0.0)) ? (i - 1) : i;
        int   const im  = clamp_idx(im0, phi_lo.x, phi_hi.x - 1);
        int   const jc  = clamp_idx(j, phi_lo.y, phi_hi.y - 1);
        amrex::Real const w = (u < amrex::Real(0.0)) ? (u + amrex::Real(1.0)) : u;
        amrex::Real const pj  = phi(im, jc,     kk) * (amrex::Real(1.0) - w)
                              + phi(im + 1, jc, kk) * w;
        amrex::Real const pj1 = phi(im, jc + 1,     kk) * (amrex::Real(1.0) - w)
                              + phi(im + 1, jc + 1, kk) * w;
        amrex::Real const pc = amrex::Real(0.5) * (pj + pj1);
        if (pc > amrex::Real(0.0)) { all_open = false; }
        if (s > 0) {
            len += open_fraction_1d(prev, pc) * h / static_cast<amrex::Real>(nseg);
        }
        prev = pc;
    }
    // exact shortcut on regular regions: guarantees bitwise h (needed by L2a)
    return all_open ? h : len;
}

/** \brief Open length along z on the dual line through (x_{i+1/2}, z_j). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real open_len_z_from_ls (amrex::Array4<amrex::Real const> const& phi,
                                amrex::Dim3 const phi_lo, amrex::Dim3 const phi_hi,
                                int i, int j, int k, amrex::Real h, int nseg)
{
    int const kk = clamp_idx(k, phi_lo.z, phi_hi.z);
    amrex::Real len = amrex::Real(0.0);
    amrex::Real prev = amrex::Real(0.0);
    bool all_open = true;
    for (int s = 0; s <= nseg; ++s) {
        amrex::Real const u = -amrex::Real(0.5)
            + static_cast<amrex::Real>(s) / static_cast<amrex::Real>(nseg);
        int   const jm0 = (u < amrex::Real(0.0)) ? (j - 1) : j;
        int   const jm  = clamp_idx(jm0, phi_lo.y, phi_hi.y - 1);
        int   const ic  = clamp_idx(i, phi_lo.x, phi_hi.x - 1);
        amrex::Real const w = (u < amrex::Real(0.0)) ? (u + amrex::Real(1.0)) : u;
        amrex::Real const pi  = phi(ic,     jm, kk) * (amrex::Real(1.0) - w)
                              + phi(ic,     jm + 1, kk) * w;
        amrex::Real const pi1 = phi(ic + 1, jm, kk) * (amrex::Real(1.0) - w)
                              + phi(ic + 1, jm + 1, kk) * w;
        amrex::Real const pc = amrex::Real(0.5) * (pi + pi1);
        if (pc > amrex::Real(0.0)) { all_open = false; }
        if (s > 0) {
            len += open_fraction_1d(prev, pc) * h / static_cast<amrex::Real>(nseg);
        }
        prev = pc;
    }
    return all_open ? h : len;
}

} // anonymous namespace


void
web::MarkReducedShapeCells (
    std::unique_ptr<amrex::iMultiFab> & eb_reduce_particle_shape,
    amrex::EBFArrayBoxFactory const & eb_fact,
    int const particle_shape_order,
    const amrex::Periodicity& periodicity)
{
    // Pre-fill array with 0, including in the ghost cells outside of the domain.
    // (The guard cells in the domain will be updated by `FillBoundary` at the end of this function.)
    eb_reduce_particle_shape->setVal(0, eb_reduce_particle_shape->nGrow());

    // Extract structures for embedded boundaries
    amrex::FabArray<amrex::EBCellFlagFab> const& eb_flag = eb_fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*eb_reduce_particle_shape); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        amrex::Array4<int> const & eb_reduce_particle_shape_arr = eb_reduce_particle_shape->array(mfi);

        // Check if the box (including one layer of guard cells) contains a mix of covered and regular cells
        const amrex::Box eb_info_box = mfi.tilebox(amrex::IntVect::TheCellVector()).grow(1);
        amrex::FabType const fab_type = eb_flag[mfi].getType( eb_info_box );

        if (fab_type == amrex::FabType::regular) { // All cells in the box are regular

            // Every cell in box is regular: do not reduce particle shape in any cell
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                eb_reduce_particle_shape_arr(i, j, k) = 0;
            });

        } else if (fab_type == amrex::FabType::covered) { // All cells in the box are covered

            // Every cell in box is fully covered: reduce particle shape
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                eb_reduce_particle_shape_arr(i, j, k) = 1;
            });

        } else { // The box contains a mix of covered and regular cells

            auto const & flag = eb_flag[mfi].array();

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                // Check if any of the neighboring cells over which the particle shape might extend
                // are either partially or fully covered. In this case, set eb_reduce_particle_shape_arr
                // to one for this cell, to indicate that the particle should use an order 1 shape
                // (This ensures that the particle never deposits any charge in a partially or
                // fully covered cell, even with higher-order shapes)
                // Note: in the code below `particle_shape_order/2` corresponds to the number of neighboring cells
                // over which the shape factor could extend, in each direction.
                int const i_start = i-particle_shape_order/2;
                int const i_end = i+particle_shape_order/2;
#if AMREX_SPACEDIM > 1
                int const j_start = j-particle_shape_order/2;
                int const j_end = j+particle_shape_order/2;
#else
                int const j_start = j;
                int const j_end = j;
#endif
#if AMREX_SPACEDIM > 2
                int const k_start = k-particle_shape_order/2;
                int const k_end = k+particle_shape_order/2;
#else
                int const k_start = k;
                int const k_end = k;
#endif
                int reduce_shape = 0;
                for (int i_cell = i_start; i_cell <= i_end; ++i_cell) {
                    for (int j_cell = j_start; j_cell <= j_end; ++j_cell) {
                        for (int k_cell = k_start; k_cell <= k_end; ++k_cell) {
                            // `isRegular` returns `false` if the cell is either partially or fully covered.
                            if ( !flag(i_cell, j_cell, k_cell).isRegular() ) {
                                reduce_shape = 1;
                            }
                        }
                    }
                }
                eb_reduce_particle_shape_arr(i, j, k) = reduce_shape;
            });

        }

    }
    // FillBoundary to set the values in the guard cells
    eb_reduce_particle_shape->FillBoundary(periodicity);
}

void
web::MarkUpdateCellsStairCase (
    std::array< std::unique_ptr<amrex::iMultiFab>,3> & eb_update,
    ablastr::fields::VectorField const& field,
    amrex::EBFArrayBoxFactory const & eb_fact,
    const amrex::Periodicity& periodicity )
{

    // Extract structures for embedded boundaries
    amrex::FabArray<amrex::EBCellFlagFab> const& eb_flag = eb_fact.getMultiEBCellFlagFab();

    for (int idim = 0; idim < 3; ++idim) {

        // Pre-fill the flags with 1 (i.e. "update this point"), including the
        // ghost cells outside of the domain: guard cells beyond a non-periodic
        // domain boundary are not covered by the valid-region marking below,
        // nor by the final `FillBoundary`, but they are read by consumers that
        // loop over grown tileboxes (e.g. `CalculateCurrentAmpere` or
        // `ComputeExternalFieldOnGridUsingParser`). Pre-filling here (rather
        // than only at allocation) also keeps the flags well-defined when they
        // are re-allocated during load balancing.
        // (The guard cells in the domain will be updated by `FillBoundary` at
        // the end of this function.)
        eb_update[idim]->setVal(1, eb_update[idim]->nGrow());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*field[idim]); mfi.isValid(); ++mfi) {

            const amrex::Box& box = mfi.tilebox();
            amrex::Array4<int> const & eb_update_arr = eb_update[idim]->array(mfi);

            // Check if the box (including one layer of guard cells) contains a mix of covered and regular cells
            const amrex::Box eb_info_box = mfi.tilebox(amrex::IntVect::TheCellVector()).grow(1);
            amrex::FabType const fab_type = eb_flag[mfi].getType( eb_info_box );

            if (fab_type == amrex::FabType::regular) { // All cells in the box are regular

                // Every cell in box is regular: update field in every cell
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    eb_update_arr(i, j, k) = 1;
                });

            } else if (fab_type == amrex::FabType::covered) { // All cells in the box are covered

                // Every cell in box is fully covered: do not update field
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    eb_update_arr(i, j, k) = 0;
                });

            } else { // The box contains a mix of covered and regular cells

                auto const & flag = eb_flag[mfi].array();
                auto index_type = field[idim]->ixType();

                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                    // Stair-case approximation: If neighboring cells of this gridpoint
                    // are either partially or fully covered: do not update field

                    // The number of cells that we need to check depend on the index type
                    // of the `eb_update_arr` in each direction.
                    // If `eb_update_arr` is nodal in a given direction, we need to check the cells
                    // to the left and right of this nodal gridpoint.
                    // For instance, if `eb_update_arr` is nodal in the first dimension, we need
                    // to check the cells at index i-1 and at index i, since, with AMReX indexing conventions,
                    // these are the neighboring cells for the nodal gripoint at index i.
                    // If `eb_update_arr` is cell-centerd in a given direction, we only need to check
                    // the cell at the same position (e.g., in the first dimension: the cell at index i).
                    int const i_start = ( index_type.nodeCentered(0) )? i-1 : i;
#if AMREX_SPACEDIM > 1
                    int const j_start = ( index_type.nodeCentered(1) )? j-1 : j;
#else
                    int const j_start = j;
#endif
#if AMREX_SPACEDIM > 2
                    int const k_start = ( index_type.nodeCentered(2) )? k-1 : k;
#else
                    int const k_start = k;
#endif
                    // Loop over neighboring cells
                    int eb_update_flag = 1;
                    for (int i_cell = i_start; i_cell <= i; ++i_cell) {
                        for (int j_cell = j_start; j_cell <= j; ++j_cell) {
                            for (int k_cell = k_start; k_cell <= k; ++k_cell) {
                                // If one of the neighboring is either partially or fully covered
                                // (i.e. if they are not regular cells), do not update field
                                // (`isRegular` returns `false` if the cell is either partially or fully covered.)
                                if ( !flag(i_cell, j_cell, k_cell).isRegular() ) {
                                    eb_update_flag = 0;
                                }
                            }
                        }
                    }
                    eb_update_arr(i, j, k) = eb_update_flag;
                });

            }

        }
        // Populate guard cells
        eb_update[idim]->FillBoundary(periodicity);
    }

}

void
web::MarkUpdateECellsECT (
    std::array< std::unique_ptr<amrex::iMultiFab>,3> & eb_update_E,
    ablastr::fields::VectorField const& edge_lengths,
    amrex::MultiFab const* dual_A_y,
    amrex::Real A_stab )
{

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(*eb_update_E[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // z is degenerate in 2D: keep the loop at k = 0, whatever ghost the
        // iMultiFab carries in 3D (otherwise the EB arrays are indexed at k < 0).
        amrex::IntVect const ngZ = eb_update_E[0]->nGrowVect();
        amrex::IntVect const ng2dZ(ngZ[0], ngZ[1]);
        const amrex::Box& tbx = mfi.tilebox( eb_update_E[0]->ixType().toIntVect(), ng2dZ );
        const amrex::Box& tby = mfi.tilebox( eb_update_E[1]->ixType().toIntVect(), ng2dZ );
        const amrex::Box& tbz = mfi.tilebox( eb_update_E[2]->ixType().toIntVect(), ng2dZ );

        amrex::Array4<int> const & eb_update_Ex_arr = eb_update_E[0]->array(mfi);
        amrex::Array4<int> const & eb_update_Ey_arr = eb_update_E[1]->array(mfi);
        amrex::Array4<int> const & eb_update_Ez_arr = eb_update_E[2]->array(mfi);

        amrex::Array4<amrex::Real> const & lx_arr = edge_lengths[0]->array(mfi);
        amrex::Array4<amrex::Real> const & lz_arr = edge_lengths[2]->array(mfi);
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_XZ)
        // In XZ, ly is the openness of the mesh node on which Ey lives.
        amrex::Array4<amrex::Real> const & ly_arr = edge_lengths[1]->array(mfi);
#endif
#if defined(WARPX_DIM_XZ)
        // When the conformal-E dual geometry is available, use the same quantity that
        // the kernel divides by (A_y) to decide whether a node is "covered". The ly
        // criterion also zeroes nodes with only one covered side, which have A_y > 0,
        // i.e. it freezes part of the cut nodes for no reason (C_SIDE_FIX_PLAN.md F3).
        bool const has_dual_Ay = (dual_A_y != nullptr);
        amrex::Array4<amrex::Real const> dual_Ay_arr;
        if (has_dual_Ay) {
            dual_Ay_arr = dual_A_y->const_array(mfi);
        }
#endif
#if defined(WARPX_DIM_RZ)
        amrex::Dim3 const lx_lo = amrex::lbound(lx_arr);
        amrex::Dim3 const lx_hi = amrex::ubound(lx_arr);
        amrex::Dim3 const lz_lo = amrex::lbound(lz_arr);
        amrex::Dim3 const lz_hi = amrex::ubound(lz_arr);
#endif

        amrex::ParallelFor (tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Do not update Ex if the edge on which it lives is fully covered
                eb_update_Ex_arr(i, j, k) = (lx_arr(i, j, k) == 0)? 0 : 1;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#if defined(WARPX_DIM_3D)
                // In 3D: Do not update Ey if the edge on which it lives is fully covered
                eb_update_Ey_arr(i, j, k) = (ly_arr(i, j, k) == 0)? 0 : 1;
#elif defined(WARPX_DIM_XZ)
                // In XZ: Ey is associated with a mesh node. Prefer the dual area A_y
                // (the kernel divides by it); fall back to ly when it is not available.
                if (has_dual_Ay) {
                    eb_update_Ey_arr(i, j, k) = (dual_Ay_arr(i, j, k) <= A_stab)? 0 : 1;
                } else {
                    eb_update_Ey_arr(i, j, k) = (ly_arr(i, j, k) == 0)? 0 : 1;
                }
#elif defined(WARPX_DIM_RZ)
                // In RZ: Ey is associated with a mesh node,
                // so we need to check if  the mesh node is covered
                if((lx_arr(std::min(i  , lx_hi.x), std::min(j  , lx_hi.y), k)==0)
                 ||(lx_arr(std::max(i-1, lx_lo.x), std::min(j  , lx_hi.y), k)==0)
                 ||(lz_arr(std::min(i  , lz_hi.x), std::min(j  , lz_hi.y), k)==0)
                 ||(lz_arr(std::min(i  , lz_hi.x), std::max(j-1, lz_lo.y), k)==0)) {
                    eb_update_Ey_arr(i, j, k) = 0;
                } else {
                    eb_update_Ey_arr(i, j, k) = 1;
                }
#endif
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Do not update Ez if the edge on which it lives is fully covered
                eb_update_Ez_arr(i, j, k) = (lz_arr(i, j, k) == 0)? 0 : 1;
            }
        );

    }
}

void
web::MarkUpdateBCellsECT (
    std::array< std::unique_ptr<amrex::iMultiFab>,3> & eb_update_B,
    ablastr::fields::VectorField const& face_areas,
    ablastr::fields::VectorField const& edge_lengths )
{

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(*eb_update_B[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // z is degenerate in 2D: keep the loop at k = 0, whatever ghost the
        // iMultiFab carries in 3D (otherwise the EB arrays are indexed at k < 0).
        amrex::IntVect const ngZ = eb_update_B[0]->nGrowVect();
        amrex::IntVect const ng2dZ(ngZ[0], ngZ[1]);
        const amrex::Box& tbx = mfi.tilebox( eb_update_B[0]->ixType().toIntVect(), ng2dZ );
        const amrex::Box& tby = mfi.tilebox( eb_update_B[1]->ixType().toIntVect(), ng2dZ );
        const amrex::Box& tbz = mfi.tilebox( eb_update_B[2]->ixType().toIntVect(), ng2dZ );

        amrex::Array4<int> const & eb_update_Bx_arr = eb_update_B[0]->array(mfi);
        amrex::Array4<int> const & eb_update_By_arr = eb_update_B[1]->array(mfi);
        amrex::Array4<int> const & eb_update_Bz_arr = eb_update_B[2]->array(mfi);

#ifdef WARPX_DIM_3D
        amrex::Array4<amrex::Real> const & Sx_arr = face_areas[0]->array(mfi);
        amrex::Array4<amrex::Real> const & Sy_arr = face_areas[1]->array(mfi);
        amrex::Array4<amrex::Real> const & Sz_arr = face_areas[2]->array(mfi);
        amrex::ignore_unused(edge_lengths);
#elif defined(WARPX_DIM_XZ)
        // In XZ, Bx and Bz live on the faces normal to x and to z, whose
        // per-unit-y open lengths are Sx and Sz.
        amrex::Array4<amrex::Real> const & Sx_arr = face_areas[0]->array(mfi);
        amrex::Array4<amrex::Real> const & Sy_arr = face_areas[1]->array(mfi);
        amrex::Array4<amrex::Real> const & Sz_arr = face_areas[2]->array(mfi);
        amrex::ignore_unused(edge_lengths);
#elif defined(WARPX_DIM_RZ)
        amrex::Array4<amrex::Real> const & Sy_arr = face_areas[1]->array(mfi);
        amrex::Array4<amrex::Real> const & lx_arr = edge_lengths[0]->array(mfi);
        amrex::Array4<amrex::Real> const & lz_arr = edge_lengths[2]->array(mfi);
#endif
        amrex::ParallelFor (tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef WARPX_DIM_3D
                // In 3D: do not update Bx if the face on which it lives is fully covered
                eb_update_Bx_arr(i, j, k) = (Sx_arr(i, j, k) == 0)? 0 : 1;
#elif defined(WARPX_DIM_XZ)
                // In XZ: do not update Bx if the face on which it lives is fully covered
                eb_update_Bx_arr(i, j, k) = (Sx_arr(i, j, k) == 0)? 0 : 1;
#elif defined(WARPX_DIM_RZ)
                //In RZ, Bx lives on a z-edge ; do not update if fully covered
                eb_update_Bx_arr(i, j, k) = (lz_arr(i, j, k) == 0)? 0 : 1;
#endif
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Do not update By if the face on which it lives is fully covered
                eb_update_By_arr(i, j, k) = (Sy_arr(i, j, k) == 0)? 0 : 1;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef WARPX_DIM_3D
                // In 3D: do not update Bz if the face on which it lives is fully covered
                eb_update_Bz_arr(i, j, k) = (Sz_arr(i, j, k) == 0)? 0 : 1;
#elif defined(WARPX_DIM_XZ)
                // In XZ: do not update Bz if the face on which it lives is fully covered
                eb_update_Bz_arr(i, j, k) = (Sz_arr(i, j, k) == 0)? 0 : 1;
#elif defined(WARPX_DIM_RZ)
                //In RZ, Bz lives on a x-edge ; do not update if fully covered
                eb_update_Bz_arr(i, j, k) = (lx_arr(i, j, k) == 0)? 0 : 1;
#endif
            }
        );

    }
}

void
web::MarkExtensionCells (
    const std::array<amrex::Real,3>& cell_size,
    std::array< std::unique_ptr<amrex::iMultiFab>, 3 > & flag_info_face,
    std::array< std::unique_ptr<amrex::iMultiFab>, 3 > & flag_ext_face,
    const ablastr::fields::VectorField& b_field,
    const ablastr::fields::VectorField& face_areas,
    const ablastr::fields::VectorField& edge_lengths,
    const ablastr::fields::VectorField& area_mod)
{
    using ablastr::fields::Direction;
    using warpx::fields::FieldType;

#ifdef WARPX_DIM_RZ
    amrex::ignore_unused(cell_size, flag_info_face, flag_ext_face, b_field,
        face_areas, edge_lengths, area_mod);

#elif !defined(WARPX_DIM_3D) && !defined(WARPX_DIM_XZ)

    WARPX_ABORT_WITH_MESSAGE("MarkExtensionCells only implemented in 2D and 3D");

#else

    for (int idim = 0; idim < 3; ++idim) {

        // In XZ all three components are classified: idim == 0 and idim == 2
        // (the faces of the TE mode) are stabilized with the BCK correction only,
        // see the XZ branch in the kernel below.
        for (amrex::MFIter mfi(*b_field[idim]); mfi.isValid(); ++mfi) {
            auto* face_areas_idim_max_lev = face_areas[idim];

            const amrex::Box& box = mfi.tilebox(face_areas_idim_max_lev->ixType().toIntVect(),
                                                face_areas_idim_max_lev->nGrowVect() );

            auto const& S = face_areas_idim_max_lev->array(mfi);
            auto const& flag_info_face_data = flag_info_face[idim]->array(mfi);
            auto const& flag_ext_face_data = flag_ext_face[idim]->array(mfi);
            auto const& lx = edge_lengths[0]->array(mfi);
            auto const& ly = edge_lengths[1]->array(mfi);
            auto const& lz = edge_lengths[2]->array(mfi);
            auto const& mod_areas_dim_data = area_mod[idim]->array(mfi);

            const amrex::Real dx = cell_size[0];
            const amrex::Real dy = cell_size[1];
            const amrex::Real dz = cell_size[2];

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Minimal area for this cell to be stable
                mod_areas_dim_data(i, j, k) = S(i, j, k);
                double S_stab;
                if (idim == 0){
#  if defined(WARPX_DIM_XZ)
                    // Per-unit-y: the face spans dz and its two end points are the
                    // nodes (i, j) and (i, j + 1); the z-edge at x = i lies on it.
                    S_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j + 1, k) * dz,
                                             lz(i, j, k) * dy});
#  else
                    S_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                    lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
#  endif
                }else if (idim == 1){

#  if defined(WARPX_DIM_XZ)
                    S_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j + 1, k) * dz,
                                            lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
#  else
                    S_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                            lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
#  endif
                }else {
#  if defined(WARPX_DIM_XZ)
                    // Per-unit-y: the face spans dx and its two end points are the
                    // nodes (i, j) and (i + 1, j); the x-edge at z = j lies on it.
                    S_stab = 0.5 * std::max({ly(i, j, k) * dx, ly(i + 1, j, k) * dx,
                                             lx(i, j, k) * dy});
#  else
                    S_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                             ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
#  endif
                }

// TE faces (idim == 0 and 2 in XZ) can borrow area from their in-plane
                // neighbours (G2). With algo.ect_te_borrowing = 0 they fall back to the
                // pre-G2 behaviour: no extension, tiny faces enlarged by BCK.
                if ((idim == 0 || idim == 2) && (WarpX::ect_te_borrowing == 0)) {
                    flag_ext_face_data(i, j, k) = 0;
                    if (S(i, j, k) > 0) {
                        if (S(i, j, k) < S_stab) {
                            S(i, j, k) = S_stab;
                            mod_areas_dim_data(i, j, k) = S_stab;
                            flag_info_face_data(i, j, k) = FaceInfo::bck_stabilized;
                        } else {
                            flag_info_face_data(i, j, k) = FaceInfo::available;
                        }
                    }
                    return;
                }

                // Does this face need to be extended?
                // The difference between flag_info_face and flag_ext_face is that:
                //     - for every face flag_info_face contains one of the FaceInfo::Flag values.
                //       Here we only take care of FaceInfo::extended and FaceInfo::available. The
                //       entries corresponding to the intruded faces are going to be set in the
                //       function ComputeFaceExtensions
                //     - for every face flag_ext_face contains a:
                //          * 1 if the face needs to be extended
                //          * 0 otherwise
                //       In the function ComputeFaceExtensions, after the cells are extended, the
                //       corresponding entries in flag_ext_face are set to zero. This helps to keep
                //       track of which cells could not be extended
                flag_ext_face_data(i, j, k) = int(S(i, j, k) < S_stab && S(i, j, k) > 0);
                if(flag_ext_face_data(i, j, k)){
                    flag_info_face_data(i, j, k) = FaceInfo::extended;
                }
                // Is this face available to lend area to other faces?
                // The criterion is that the face has to be interior and not already unstable itself
                if(int(S(i, j, k) > 0 && !flag_ext_face_data(i, j, k))) {
                    flag_info_face_data(i, j, k) = FaceInfo::available;
                }
            });
        }
    }
#endif
}

void
web::ComputeEdgeLengths (
    ablastr::fields::VectorField& edge_lengths,
    const amrex::EBFArrayBoxFactory& eb_fact)
{
    BL_PROFILE("ComputeEdgeLengths");

#if !defined(WARPX_DIM_3D) && !defined(WARPX_DIM_XZ) && !defined(WARPX_DIM_RZ)
    WARPX_ABORT_WITH_MESSAGE("ComputeEdgeLengths only implemented in 2D and 3D");
#endif

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &edge_centroid = eb_fact.getEdgeCent();
    for (int idim = 0; idim < 3; ++idim){
#if defined(WARPX_DIM_XZ)
        // In XZ the y-directed edge degenerates into a mesh node, which has no
        // edge centroid: it is filled in after lx and lz have been computed.
        if (idim == 1) { continue; }
#elif defined(WARPX_DIM_RZ)
        if (idim == 1) {
            edge_lengths[1]->setVal(0.);
            continue;
        }
#endif
        for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi){
            amrex::Box const box = mfi.tilebox(edge_lengths[idim]->ixType().toIntVect(),
                                               edge_lengths[idim]->nGrowVect());
            amrex::FabType const fab_type = flags[mfi].getType(box);
            auto const &edge_lengths_dim = edge_lengths[idim]->array(mfi);

            if (fab_type == amrex::FabType::regular) {
                // every cell in box is all regular
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 1.;
                });
            } else if (fab_type == amrex::FabType::covered) {
                // every cell in box is all covered
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 0.;
                });
            } else {
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                int idim_amrex = idim;
                if (idim == 2) { idim_amrex = 1; }
                auto const &edge_cent = edge_centroid[idim_amrex]->const_array(mfi);
#elif defined(WARPX_DIM_3D)
                auto const &edge_cent = edge_centroid[idim]->const_array(mfi);
#endif
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (edge_cent(i, j, k) == amrex::Real(-1.0)) {
                        // This edge is all covered
                        edge_lengths_dim(i, j, k) = 0.;
                    } else if (edge_cent(i, j, k) == amrex::Real(1.0)) {
                        // This edge is all open
                        edge_lengths_dim(i, j, k) = 1.;
                    } else {
                        // This edge is cut.
                        edge_lengths_dim(i, j, k) = 1 - amrex::Math::abs(amrex::Real(2.0)
                                                                        * edge_cent(i, j, k));
                    }

                });
            }
        }
    }

#if defined(WARPX_DIM_XZ)
    // Per-unit-y reduction: the y-directed edge degenerates into a mesh node, so
    // its "length" is the openness of the node itself (0 or 1). A node is
    // considered closed as soon as one of the four edges it belongs to is fully
    // covered. This mirrors the criterion used in MarkUpdateECellsECT.
    for (amrex::MFIter mfi(*edge_lengths[1]); mfi.isValid(); ++mfi){
        amrex::Box const box = mfi.tilebox(edge_lengths[1]->ixType().toIntVect(),
                                           edge_lengths[1]->nGrowVect());
        auto const &ly = edge_lengths[1]->array(mfi);
        auto const &lx = edge_lengths[0]->const_array(mfi);
        auto const &lz = edge_lengths[2]->const_array(mfi);
        amrex::Dim3 const lx_lo = amrex::lbound(lx);
        amrex::Dim3 const lx_hi = amrex::ubound(lx);
        amrex::Dim3 const lz_lo = amrex::lbound(lz);
        amrex::Dim3 const lz_hi = amrex::ubound(lz);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            bool const closed =
                   (lx(std::min(i  , lx_hi.x), std::min(j  , lx_hi.y), k) == 0.)
                || (lx(std::max(i-1, lx_lo.x), std::min(j  , lx_hi.y), k) == 0.)
                || (lz(std::min(i  , lz_hi.x), std::min(j  , lz_hi.y), k) == 0.)
                || (lz(std::min(i  , lz_hi.x), std::max(j-1, lz_lo.y), k) == 0.);
            ly(i, j, k) = closed ? amrex::Real(0.) : amrex::Real(1.);
        });
    }
#endif
}


void
web::ComputeFaceAreas (
    ablastr::fields::VectorField& face_areas,
    const amrex::EBFArrayBoxFactory& eb_fact)
{
    BL_PROFILE("ComputeFaceAreas");

#if !defined(WARPX_DIM_3D) && !defined(WARPX_DIM_XZ) && !defined(WARPX_DIM_RZ)
    WARPX_ABORT_WITH_MESSAGE("ComputeFaceAreas only implemented in 2D and 3D");
#endif

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
#if defined(WARPX_DIM_XZ)
    //In 2D the volume frac is actually the area frac. It is the openness of the
    //face normal to y, which in XZ is a whole cell face (per unit length in y).
    auto const &area_frac = eb_fact.getVolFrac();
    // Openness of the faces normal to x (entry 0) and to z (entry 1). In XZ the
    // second AMReX dimension is the physical z; both are per unit length in y.
    auto const &area_frac_nd = eb_fact.getAreaFrac();
#elif defined(WARPX_DIM_RZ)
    //In 2D the volume frac is actually the area frac.
    auto const &area_frac = eb_fact.getVolFrac();
#elif defined(WARPX_DIM_3D)
    auto const &area_frac = eb_fact.getAreaFrac();
#endif

    for (int idim = 0; idim < 3; ++idim) {
#if defined(WARPX_DIM_RZ)
        if (idim == 0 || idim == 2) {
            face_areas[idim]->setVal(0.);
            continue;
        }
#endif
#if defined(WARPX_DIM_XZ)
        if (idim == 0 || idim == 2) {
            // Per-unit-y reduction: the face normal to x (idim == 0) or to z
            // (idim == 2) is a line segment, whose openness is directly provided
            // by the AMReX area fractions. Note that a MultiCutFab has no data on
            // regular and covered boxes, so those have to be handled separately.
            int const idim_amrex = (idim == 0) ? 0 : 1;
            auto const * areafrac_nd = area_frac_nd[idim_amrex];
            for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
                amrex::Box const box = mfi.tilebox(face_areas[idim]->ixType().toIntVect(),
                                                   face_areas[idim]->nGrowVect());
                amrex::FabType const fab_type = flags[mfi].getType(box);
                auto const &face_areas_dim = face_areas[idim]->array(mfi);
                if (fab_type == amrex::FabType::regular) {
                    // every cell in box is all regular
                    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        face_areas_dim(i, j, k) = amrex::Real(1.);
                    });
                } else if (fab_type == amrex::FabType::covered) {
                    // every cell in box is all covered
                    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        face_areas_dim(i, j, k) = amrex::Real(0.);
                    });
                } else {
                    auto const &area_frac_dim = areafrac_nd->const_array(mfi);
                    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        face_areas_dim(i, j, k) = area_frac_dim(i, j, k);
                    });
                }
            }
            continue;
        }
#endif
        for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
            amrex::Box const box = mfi.tilebox(face_areas[idim]->ixType().toIntVect(),
                                               face_areas[idim]->nGrowVect());
            amrex::FabType const fab_type = flags[mfi].getType(box);
            auto const &face_areas_dim = face_areas[idim]->array(mfi);
            if (fab_type == amrex::FabType::regular) {
                // every cell in box is all regular
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(1.);
                });
            } else if (fab_type == amrex::FabType::covered) {
                // every cell in box is all covered
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(0.);
                });
            } else {
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                auto const &face = area_frac.const_array(mfi);
#elif defined(WARPX_DIM_3D)
                auto const &face = area_frac[idim]->const_array(mfi);
#endif
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    face_areas_dim(i, j, k) = face(i, j, k);
                });
            }
        }
    }
}

void
web::ScaleEdges (
    ablastr::fields::VectorField& edge_lengths,
    const std::array<amrex::Real,3>& cell_size)
{
    BL_PROFILE("ScaleEdges");

#if !defined(WARPX_DIM_3D) && !defined(WARPX_DIM_XZ) && !defined(WARPX_DIM_RZ)
    WARPX_ABORT_WITH_MESSAGE("ScaleEdges only implemented in 2D and 3D");
#endif

    for (int idim = 0; idim < 3; ++idim){
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        if (idim == 1) { continue; }
#endif
        for (amrex::MFIter mfi(*edge_lengths[0]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.tilebox(edge_lengths[idim]->ixType().toIntVect(),
                                                edge_lengths[idim]->nGrowVect() );
            auto const &edge_lengths_dim = edge_lengths[idim]->array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edge_lengths_dim(i, j, k) *= cell_size[idim];
            });
        }
    }
}


void
web::ScaleAreas (
    ablastr::fields::VectorField& face_areas,
    const std::array<amrex::Real,3>& cell_size)
{
    BL_PROFILE("ScaleAreas");

#if !defined(WARPX_DIM_3D) && !defined(WARPX_DIM_XZ) && !defined(WARPX_DIM_RZ)
    WARPX_ABORT_WITH_MESSAGE("ScaleAreas only implemented in 2D and 3D");
#endif

    for (int idim = 0; idim < 3; ++idim) {
#if defined(WARPX_DIM_RZ)
        if (idim == 0 || idim == 2) { continue; }
#endif
        // In XZ, cell_size[1] == 1, so that full_area is the open length of the
        // face per unit length in y (i.e. dz for idim == 0 and dx for idim == 2).
        for (amrex::MFIter mfi(*face_areas[0]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.tilebox(face_areas[idim]->ixType().toIntVect(),
                                                face_areas[idim]->nGrowVect() );
            amrex::Real const full_area = cell_size[(idim+1)%3]*cell_size[(idim+2)%3];
            auto const &face_areas_dim = face_areas[idim]->array(mfi);

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                face_areas_dim(i, j, k) *= full_area;
            });

        }
    }
}


void
web::ComputeDualGeometry (
    amrex::MultiFab& dual_A_y,
    amrex::MultiFab& dual_L_x,
    amrex::MultiFab& dual_L_z,
    const amrex::EBFArrayBoxFactory& eb_fact,
    const std::array<amrex::Real,3>& cell_size,
    int scheme,
    amrex::Real A_stab)
{
    BL_PROFILE("ComputeDualGeometry");
    amrex::Print() << "[DIAG] ComputeDualGeometry called, scheme=" << scheme << "\n";

    amrex::ignore_unused(dual_A_y, dual_L_x, dual_L_z, eb_fact, cell_size, scheme, A_stab);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        scheme == 1 || scheme == 2 || scheme == 4 || scheme == 6 || scheme == 7 || scheme == 8 || scheme == 9,
        "ComputeDualGeometry: only scheme 1 (constant), 2 (volfrac) and 4 (level set) "
        "are implemented");

#if defined(WARPX_DIM_XZ)
    amrex::Real const dx = cell_size[0];
    amrex::Real const dz = cell_size[2];

    auto const& flags = eb_fact.getMultiEBCellFlagFab();
    auto const& volfrac = eb_fact.getVolFrac();

    // NOTE on indexing (see Docs/research_notes/C_SIDE_DESIGN2.md):
    //   dual_A_y is centred on the nodes  (IntVect = (1,1)) -- same place as Ey
    //   dual_L_x is centred where Bx lives (IntVect = (1,0)) -- open length along x
    //   dual_L_z is centred where Bz lives (IntVect = (0,1)) -- open length along z
    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {

        // The dual arrays only exist at k = 0 in 2D, but ng_FieldSolver carries ghost
        // cells in all three directions. Using nGrowVect() directly would make the
        // write box span k = -2..2 and index the (k=0 only) multiFabs out of bounds.
        amrex::IntVect const ngA = dual_A_y.nGrowVect();
        amrex::IntVect const ng2d(ngA[0], ngA[1]);
        amrex::Box const& bAy = mfi.tilebox(dual_A_y.ixType().toIntVect(), ng2d);
        amrex::Box const& bLx = mfi.tilebox(dual_L_x.ixType().toIntVect(), ng2d);
        amrex::Box const& bLz = mfi.tilebox(dual_L_z.ixType().toIntVect(), ng2d);

        auto const& Ay = dual_A_y.array(mfi);
        auto const& Lx = dual_L_x.array(mfi);
        auto const& Lz = dual_L_z.array(mfi);

        amrex::FabType const fab_type =
            (scheme == 1) ? amrex::FabType::regular : flags[mfi].getType(mfi.validbox());

        if (fab_type == amrex::FabType::regular) {

            amrex::ParallelFor(bAy, bLx, bLz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) { Ay(i,j,k) = dx*dz; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) { Lx(i,j,k) = dx; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) { Lz(i,j,k) = dz; });

        } else if (fab_type == amrex::FabType::covered) {

            amrex::ParallelFor(bAy, bLx, bLz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) { Ay(i,j,k) = amrex::Real(0.0); },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) { Lx(i,j,k) = amrex::Real(0.0); },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) { Lz(i,j,k) = amrex::Real(0.0); });

        } else {

            auto const& v = volfrac.const_array(mfi);

            // dual_A_y: four quadrants, one quarter of each cell.
            // A_y is a *denominator* in the kernel, so a cut node with A_y -> 0 would
            // amplify Ey by 1/A_y (hot spots inside the PEC). Raise (0, A_stab) to
            // A_stab, but keep exactly zero (fully covered) at zero so that the mask
            // still freezes it.
            // The write box includes ghost cells, so `i-1`/`j-1` can fall outside the
            // volume-fraction MultiFab (it carries fewer ghost layers). Clamp the
            // indices: this only affects the ghost region, never the valid cells.
            amrex::Dim3 const v_lo = amrex::lbound(v);
            amrex::Dim3 const v_hi = amrex::ubound(v);
            amrex::ParallelFor(bAy, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                int const im = (i > v_lo.x) ? (i-1) : v_lo.x;
                int const jm = (j > v_lo.y) ? (j-1) : v_lo.y;
                int const ip = (i < v_hi.x) ? i : v_hi.x;
                int const jp = (j < v_hi.y) ? j : v_hi.y;
                amrex::Real a = (scheme == 7) ? (dx*dz) :
                    amrex::Real(0.25)*dx*dz*
                    ( v(im,jm,k) + v(ip,jm,k) + v(im,jp,k) + v(ip,jp,k) );
                if (a > amrex::Real(0.0) && a < A_stab) { a = A_stab; }
                Ay(i,j,k) = a;
            });

            if (scheme == 4) {

                // scheme C: integrate the nodal level set along the dual lines.
                // Dimension-consistent (a 1D integral of a 1D quantity), unlike the
                // volume-fraction average of scheme B; expected second order.
                constexpr int nseg = 8;
                auto const& phi = eb_fact.getLevelSet().const_array(mfi);
                amrex::Dim3 const phi_lo = amrex::lbound(phi);
                amrex::Dim3 const phi_hi = amrex::ubound(phi);

                amrex::ParallelFor(bLx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Lx(i,j,k) = open_len_x_from_ls(phi, phi_lo, phi_hi, i, j, k, dx, nseg);
                });

                amrex::ParallelFor(bLz, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Lz(i,j,k) = open_len_z_from_ls(phi, phi_lo, phi_hi, i, j, k, dz, nseg);
                });

            } else {

                // dual_L_x(i,j): open length along x on the line z = z_{j+1/2}
                amrex::ParallelFor(bLx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    int const imx = (i > v_lo.x) ? (i-1) : v_lo.x;
                    int const ipx = (i < v_hi.x) ? i : v_hi.x;
                    amrex::Real const lx_v = amrex::Real(0.5)*dx*( v(imx,j,k) + v(ipx,j,k) );
                    amrex::Real const lx_vmin = (v(imx,j,k) < v(ipx,j,k)) ? v(imx,j,k) : v(ipx,j,k);
                    Lx(i,j,k) = (scheme == 6) ? dx
                              : (scheme == 8) ? ((lx_v > amrex::Real(0.5)*dx) ? lx_v : amrex::Real(0.5)*dx)
                              : ((scheme == 2 || scheme == 9) && lx_vmin < amrex::Real(0.2)) ? dx
                              : lx_v;
                });

                // dual_L_z(i,j): open length along z on the line x = x_{i+1/2}
                amrex::ParallelFor(bLz, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    int const jmy = (j > v_lo.y) ? (j-1) : v_lo.y;
                    int const jpy = (j < v_hi.y) ? j : v_hi.y;
                    int const ipy = (i < v_hi.x) ? i : v_hi.x;
                    amrex::Real const lz_v = amrex::Real(0.5)*dz*( v(ipy,jmy,k) + v(ipy,jpy,k) );
                    amrex::Real const lz_vmin = (v(ipy,jmy,k) < v(ipy,jpy,k)) ? v(ipy,jmy,k) : v(ipy,jpy,k);
                    Lz(i,j,k) = (scheme == 6) ? dz
                              : (scheme == 8) ? ((lz_v > amrex::Real(0.5)*dz) ? lz_v : amrex::Real(0.5)*dz)
                              : ((scheme == 2 || scheme == 9) && lz_vmin < amrex::Real(0.2)) ? dz
                              : lz_v;
                });
            }
        }
    }
#endif
    amrex::Print() << "[DIAG] dual geom: A_y=" << dual_A_y.min(0) << ".." << dual_A_y.max(0) << " Lx=" << dual_L_x.min(0) << ".." << dual_L_x.max(0) << " Lz=" << dual_L_z.min(0) << ".." << dual_L_z.max(0) << "\n";
}

#endif
