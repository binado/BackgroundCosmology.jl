# BackgroundCosmology.jl

[![CI](https://github.com/binado/BackgroundCosmology.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/binado/BackgroundCosmology.jl/actions/workflows/CI.yml)

Flat FLRW background cosmology for gravitational-wave population work, with a
[Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl)-style API and **no units**:
`H0` is in km/s/Mpc, distances in Mpc, times in Gyr. Supports flat ΛCDM, wCDM and w0waCDM
(no radiation, no curvature).

All integrals use fixed-order Gauss–Legendre rules (FastGaussQuadrature), never adaptive
quadrature, so [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) duals propagate
through every call. Results match the Python `astrogwb` stack to rounding error.

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
luminosity_dist(c, 1.0)          # Mpc
comoving_radial_dist(c, 1.0)     # Mpc
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
