# C03 — medium-solver test suite (quantitative)

A ladder from **math unit tests** (formulas, no solver) to **2D simulation experiments**
(field/interface). Every entry names its reference (analytic / independent Python /
self-consistency), its acceptance criterion, and its cost.

| # | case | feature | reference | criterion | cost |
|---|------|---------|-----------|-----------|------|
| M1 | `math/test_reference_math.py` | isotropic `chi^-1` | harmonic/arithmetic limits | ≤1e-14 | <1 s |
| M2 | same | anisotropic `tau` (diagonal) | same limits | ≤1e-14 | <1 s |
| M3 | same | general 3D `frame` | orthonormal, e_z=n | ≤1e-14 | <1 s |
| M4 | same | `expm` M1/M2 (tensor eps, tensor sigma) | `scipy.linalg.expm` | ≤1e-13 | < 2 s |
| M5 | same | `(1-e^-x)/x` series branch | `mpmath` 50 digits | ≤1e-15 | < 1 s |
| S1 | `c02_0_pec/` (existing) | sigma relaxation | `exp(-sigma t/eps)` | ≤1e-12 | 1 min |
| S2 | `c02_0_pec/` (existing) | Gauss consistency (`rho_medium`) | internal machine precision | 4e-25 (edge 7e-13) | 1 min |
| S3 | `c02_0_pec/` (existing) | sigma=0 bitwise regression | previous run | bitwise | 2 min |
| S4 | `c02_0_pec/` (existing) | sigma -> infinity, 1 step | E -> 0 | 0, finite | < 1 min |
| S5 | `lossy_slab/` | lossy plane-wave attenuation | `exp(-Im k d)`, phase `-Re k d` | amp ≤2%, phase ≤3 deg (32 ppw) | 20-60 min |
| E1 | `c02_0_pec/verify_fresnel_rt.py` | Fresnel normal s/p | analytic | R/T ≤2% (32 ppw) | 3 min |
| E2 | `c02_0_pec/verify_oblique_rt.py` | Fresnel oblique s/p | analytic | R/T ≤2-3% (32 ppw) | 5 min |
| E3 | `tensor_interface/` | tensor eps interface, s-pol | n2 = sqrt(eps_yy) Fresnel | R/T ≤2% (32 ppw) | 5 min |
| E4 | `anisotropic_slab/` | bulk tensor, phase speeds | c/sqrt(eps_xx), c/sqrt(eps_yy) | <=0.5% | 2 min |
| E5 | `c02_0_pec/` (existing) | Brewster angle | analytic | qualitative (floor) | 5 min |
| E6 | `lossy_slab/` (slab variant) | conductor slab attenuation | `exp(-2 alpha d)` | ≤3% | 5 min |
| C1 | driver (optional) | convergence 16/32/64 ppw | O(h^2) | slope ~2 | 1-2 h |

## Usage

Math layer (no WarpX build needed):

```bash
python3 c03_medium/math/test_reference_math.py
```

2D cases: each directory has an `inputs_2d` template and an `analyze.py` that reads the
`FieldProbe` dumps (`diags/reducedfiles/<probe>.txt`, columns: step, time, x, y, z,
Ex, Ey, Ez, Bx, By, Bz, S).  Run WarpX with `WARPX_BIN=...` in the case directory, then:

```bash
python3 analyze.py diags/reducedfiles
```

## Conventions

- `my_constants.eps0 = 8.8541878128e-12` (WarpX value; do not mix CODATA).
- All acceptance numbers assume **32 ppw in the medium** unless stated.
- Lossy theory uses the exact complex wavenumber
  `k = omega sqrt(mu eps) sqrt(1 - i sigma/(omega eps))`; the phase convention is
  `arg[E(z2)/E(z1)] = -Re(k) dz` for a +z travelling `cos(omega t - k z)` wave.
