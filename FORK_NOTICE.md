# Fork notice

This repository is a fork of [WarpX](https://github.com/BLAST-WarpX/warpX).
WarpX is distributed under the BSD-3-Clause license; the upstream license and
copyright notices are retained in `LICENSE.txt`, `NOTICE.txt` and `LEGAL.txt`.

## What this fork adds

- **Tensor-medium field solver** (`em_solver_medium = macroscopic`): a Yee path
  with 3x3 permittivity, permeability and conductivity tensors, including a
  per-E-component `chi^-1` row representation and the full matrix-vector update
  `E_c = alpha E_c + beta0 * sum_d (chi^-1)_{c,d} D_d`.
- **Anisotropic subpixel smoothing** (Kottke 2008, Phys. Rev. E 77, 036611) for
  dielectric interfaces, with the general `tau` transform for anisotropic
  materials, a bulk mode for uniform anisotropic media, and the interface-frame
  inverse.
- **Matrix-exponential conductivity integrator**: `E^{n+1} = M1 E^n + M2 (curlH - J)`
  with `M1 = exp(-A dt)`, `M2 = (1-exp(-A dt)) A^-1 eps^-1`, `A = eps^-1 sigma`;
  `sigma = 0` cells keep the original Yee path bit-for-bit.
- **Validation suite and research notes** under `benchmarks/cases/c02_0_pec/` and
  `Docs/research_notes/`.

## License

- **All code and documentation in this repository: BSD-3-Clause.**
  Upstream WarpX is BSD-3-Clause; this fork's additions (tensor-medium solver,
  test harnesses, research notes) are released under the same license.
- The canonical text is in `LICENSE` (a standard BSD-3-Clause notice so that
  GitHub recognises it). The upstream notices in `LICENSE.txt`, `NOTICE.txt`
  and `LEGAL.txt` are retained unmodified.
- New files contributed by this fork carry an `SPDX-License-Identifier: BSD-3-Clause`
  header. Upstream files are left untouched.


## Copyright

- Upstream WarpX: Copyright (c) 2018, The Regents of the University of
  California, through Lawrence Berkeley National Laboratory, and Lawrence
  Livermore National Security, LLC (see `LICENSE.txt`, `NOTICE.txt`, `LEGAL.txt`).
- Contributions of this fork: Copyright (c) 2026, DLD222.
  Released under the BSD-3-Clause license (see `LICENSE`).

## Research-notes index

- `Docs/research_notes/技术路径.md` — technical path (phases P0-P5)
- `Docs/research_notes/P0_QUANTITATIVE_TESTS.md` — quantitative test plan
- `Docs/research_notes/P0_TESTS_EXECUTION_LOG.md` — execution results
- `Docs/research_notes/P1_BUG_HUNT_PROTOCOL.md` — P1 tensor bug-hunt protocol
- `Docs/research_notes/MEDIUM_CODE_AUDIT.md` — medium-solver code audit
- `Docs/research_notes/CUTCELL_GEOMETRY_*.md` — cut-cell / EB geometry notes

## Attribution

Upstream WarpX: <https://github.com/BLAST-WarpX/warpX>, BSD-3-Clause,
Copyright (c) 2018, The Regents of the University of California, through
Lawrence Berkeley National Laboratory, and Lawrence Livermore National
Security, LLC.
