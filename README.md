# a_changed_simulation_project

**An unofficial fork of [WarpX](https://github.com/BLAST-WarpX/warpX), distributed
under the BSD-3-Clause license.**

This repository is derived from WarpX. It is **not** the official WarpX
repository, and it is **not** affiliated with, endorsed by, sponsored by, or
approved by the WarpX developers, the BLAST-WarpX organisation, Lawrence
Berkeley National Laboratory, Lawrence Livermore National Security, or the
U.S. Department of Energy. No such endorsement is claimed or implied. The names
of the upstream authors and contributors are used only to preserve the required
copyright notices; they are not used to promote this fork.

- Upstream project: <https://github.com/BLAST-WarpX/warpX>
- Upstream documentation: <https://warpx.readthedocs.io>
- Upstream license: **BSD-3-Clause-LBNL** (see `LICENSE.txt`, `NOTICE.txt`, `LEGAL.txt`,
  kept unmodified in this repository)
- License of this fork: **BSD-3-Clause** (see `LICENSE`)

## Differences from upstream WarpX

**Project level.** This is an **unofficial fork** of
[WarpX](https://github.com/BLAST-WarpX/warpX). It is not the official WarpX
repository and is not affiliated with, endorsed by, sponsored by or approved
by the WarpX developers, BLAST-WarpX, LBNL, LLNL or the U.S. DOE. Upstream code
and notices are retained unmodified (`LICENSE.txt`, `NOTICE.txt`, `LEGAL.txt`);
this fork's additions are released under BSD-3-Clause (`LICENSE`).

**Technical differences** (all on the macroscopic medium path, Yee solver):

| Area | Upstream WarpX | This fork |
|---|---|---|
| Material tensors | scalar (or diagonal) epsilon / mu / sigma | 3x3 **tensor** epsilon, mu, sigma on the Yee path; per-E-component `chi^-1` rows; full `E_c = alpha E_c + beta0 sum_d (chi^-1)_{c,d} D_d` |
| Tensor input | -- | `epsilon_tensor_xx..yz`, `epsilon2_tensor_*` (6 components each); **bulk mode** for uniform anisotropic media (no interface) |
| Dielectric interface | scalar subpixel (Kottke, isotropic) | **anisotropic subpixel**: Kottke tau transform for anisotropic materials, interface-frame inverse, general 3D frame (`n_y != 0`), f/n from the interface parser |
| Conductivity | local backward-Euler / Lax-Wendroff update | **matrix-exponential integrator** `E = M1 E + M2 (curlH - J)` with `M1 = exp(-A dt)`, `M2 = (1-exp(-A dt)) A^-1 eps^-1`; `sigma = 0` cells keep the original Yee path bit-for-bit |
| Validation | upstream test suite | added quantitative suite: Fresnel (normal/oblique, s/p), lossy-wave attenuation, disk resonances, tensor unit tests (`benchmarks/cases/c02_0_pec/`), and 41 research notes (`Docs/research_notes/`) |

**What is not changed.** The EB/cut-cell (`ect`) solver keeps its scalar
`inv_epsr` path; the particle/PIC side, PSATD, and the upstream test suite are
untouched. The tensor/subpixel path is **Yee-only** (`macroscopic_subpixel` with
`maxwell_solver = ect` aborts). Known upstream limitation: PML combined with
`em_solver_medium = macroscopic` is not supported (see `Docs/research_notes/MEDIUM_CODE_AUDIT.md`).
A tensor conductivity input (`sigma_tensor_*`) is not implemented yet.

See [FORK_NOTICE.md](FORK_NOTICE.md) for scope, licensing and the research-notes index.

## License

BSD-3-Clause. The canonical text is in [LICENSE](LICENSE). The upstream notices
in `LICENSE.txt`, `NOTICE.txt` and `LEGAL.txt` are retained unmodified, as required
by the upstream license. Files added by this fork carry an
`SPDX-License-Identifier: BSD-3-Clause` header.

---

## Upstream README (reproduced below for reference)

> The content from here on is the WarpX README as inherited. Its badges and links
> point at the **upstream** project (Azure Pipelines, ReadTheDocs, CDash, Zenodo);
> they do not refer to this fork and do not indicate any endorsement of it.
### WarpX (upstream project)


## Overview

WarpX is an advanced **electromagnetic & electrostatic Particle-In-Cell** code.
It supports many features including Perfectly-Matched Layers (PML), mesh refinement, and the boosted-frame technique.

WarpX is a *highly-parallel and highly-optimized code*, which can run on GPUs and multi-core CPUs, and includes load balancing capabilities.
WarpX scales to the world's largest supercomputers and was awarded the [2022 ACM Gordon Bell Prize](https://www.exascaleproject.org/ecp-supported-collaborative-teams-win-the-2022-acm-gordon-bell-prize-and-special-prize/).

## Documentation


In order to learn how to install and run the code, please see the online documentation:
https://warpx.readthedocs.io

To contact the developers, feel free to open an issue on this repo, or visit our discussions page at https://github.com/BLAST-WarpX/warpx/discussions

## Contributing


Our workflow is described in [CONTRIBUTING.md](CONTRIBUTING.md).
We invite you to contribute to WarpX in any form following our [Code of Conduct](https://warpx.readthedocs.io/en/latest/coc.html), e.g., contribute to [discussions](https://github.com/BLAST-WarpX/warpx/discussions), help each other in [issues](https://github.com/BLAST-WarpX/warpx/issues), fix bugs, or add [documentation](https://warpx.readthedocs.io/en/latest/developers/documentation.html) and new functionality!

## Governance

WarpX is hosted by the High Performance Software Foundation (HPSF).
If your organization wants to help steer the evolution of the HPC software ecosystem, visit [hpsf.io](https://hpsf.io) and consider joining!

The WarpX open governance model is described in [GOVERNANCE.md](GOVERNANCE.md).

## Copyright Notice

WarpX Copyright (c) 2018, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Innovation & Partnerships Office at
IPO@lbl.gov.

Please see the full license agreement in [LICENSE.txt](LICENSE.txt).
Please see the notices in [NOTICE.txt](NOTICE.txt).
The SPDX license identifier is `BSD-3-Clause-LBNL`.
