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

## What this fork adds

- A tensor-medium field solver (`algo.em_solver_medium = macroscopic`): Yee path with
  3x3 permittivity / permeability / conductivity tensors, per-E-component
  `chi^-1` rows, and the full matrix-vector update `E = alpha E + beta0 sum_d chi^-1_cd D_d`.
- Kottke subpixel interface averaging (tau transform) with 6-component tensor input.
- A matrix-exponential integrator for conductive media, precomputing
  `M1 = exp(-A dt)` and `M2 = (1-exp(-A dt)) A^-1 eps^-1` per cell.
- A quantitative validation suite (`benchmarks/`, `Docs/research_notes/`).

See [FORK_NOTICE.md](FORK_NOTICE.md) for scope, the licence statement and the
research-notes index.

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
