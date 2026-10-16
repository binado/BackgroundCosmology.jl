# BackgroundCosmology.jl

[![CI](https://github.com/binado/BackgroundCosmology.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/binado/BackgroundCosmology.jl/actions/workflows/CI.yml)

Flat FLRW background cosmology for gravitational-wave population work, with a
[Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl)-style API and **no units**:
`H0` is in km/s/Mpc, distances in Mpc, times in Gyr. Supports flat ΛCDM, wCDM and w0waCDM
(no radiation, no curvature).

By default all integrals use fixed-order Gauss–Legendre rules (FastGaussQuadrature), never
adaptive quadrature, so [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) duals
propagate through every call; see [Choosing the quadrature method](#choosing-the-quadrature-method)
to change that. Results match the Python `astrogwb` stack to rounding error.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/binado/BackgroundCosmology.jl")
```

## Example

```julia
using BackgroundCosmology

c = cosmology(; H0 = 67.0, Ωm = 0.315)             # FlatLCDM; Planck 2018 if no keywords
c = cosmology(; H0 = 67.0, Ωm = 0.315, w0 = -0.9)  # FlatWCDM
c = cosmology(; H0 = 67.0, Ωm = 0.315, w0 = -0.9, wa = 0.2)  # FlatW0WaCDM

E(c, 1.0)                        # H(z) / H₀
luminosity_distance(c, 1.0)          # Mpc
comoving_radial_distance(c, 1.0)     # Mpc
comoving_volume_element(c, 1.0)  # Mpc³ per steradian per unit redshift
age(c), lookback_time(c, 1.0)    # Gyr

# All-sky distances and volume on your own grid in a single pass
z = collect(LinRange(0.0, 2.0, 101))
grid = distance_and_volume_grid(c, z)
# -> (; comoving_distance, luminosity_distance, differential_comoving_volume)
```

`distance_and_volume_grid` returns the **all-sky** `dV_c/dz` (`4π × comoving_volume_element`,
Mpc³ per unit redshift), which is what redshift-population codes such as
[GWDistributions.jl](https://github.com/binado/GWDistributions.jl) expect.

For sampler hyperparameters, build from a `NamedTuple` (extra keys are ignored):

```julia
cosmology((H0 = 67.0, Ωm = 0.315, Ξ₀ = 1.2))   # FlatLCDM(67.0, 0.315)
```

## Development

```
just install-runic   # once
just fmt             # format with Runic
just test            # Pkg.test() on Julia 1.12
```

## Choosing the quadrature method

Every integrating function (`comoving_radial_distance` and everything built on it,
`distance_and_volume_grid`, `age`, `lookback_time`) takes a `method` keyword that is
handed to `integrate(f, lo, hi, method)`. The change of variables and the grid structure
stay inside the package; `method` only integrates a smooth 1-D function.

```julia
using BackgroundCosmology
using BackgroundCosmology: GaussLegendreRule

luminosity_distance(c, 2.0; method = GaussLegendreRule(80))  # higher-order fixed rule

using Integrals                                           # activates the extension
luminosity_distance(c, 2.0; method = QuadGKJL())              # any SciML 1-D algorithm
luminosity_distance(c, 2.0; method = (QuadGKJL(), (; reltol = 1e-12)))  # with tolerances
age(c; method = Integrals.GaussLegendre(; n = 100))

# your own rule: add a method for your type
struct MyRule end
BackgroundCosmology.integrate(f, lo, hi, ::MyRule) = ...
```

The defaults are fixed rules (40 nodes for `d_c`, 64 for `age`, 4 per interval on a grid)
that are accurate to ~1e-12. Adaptive algorithms from Integrals.jl use their own, looser
default tolerances, so pass `reltol` if you need parity. Integrals.jl is an optional
dependency; `Integrals.GaussLegendre` runs on the same kernel with no overhead.
