# Repository Guidelines

BackgroundCosmology.jl: flat FLRW background cosmology (ΛCDM, wCDM, w0waCDM) with a
[Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl)-style API and **no units**:
`H0` in km/s/Mpc, distances in Mpc, times in Gyr.

- Layout: `src/model.jl` (types, `cosmology` builders, `E`/`H`/`w`), `src/quadrature.jl`
  (fixed Gauss–Legendre rules from FastGaussQuadrature), `src/distance.jl`, `src/time.jl`
  (`age`, `lookback_time`), `src/conversion.jl`.
- All quadrature is fixed-order Gauss–Legendre, never adaptive, so ForwardDiff duals
  propagate through every call. Keep new integrals that way; pin them against adaptive
  `QuadGK` (a test-only dependency) in the tests.
- `distance_and_volume_grid` returns the **all-sky** volume (`4π × comoving_volume_element`)
  for parity with the Python `astrogwb` stack; `comoving_volume_element` is per steradian.
- `cosmology(::Type, ::NamedTuple)` / `cosmology(::NamedTuple)` build from sampler
  hyperparameters and ignore extra keys; keep them.
- Format with Runic (`just install-runic` once, then `just fmt`; CI runs `just fmt-check`).
- Test with `just test` (Julia 1.12; includes Aqua, and JET on 1.12 only).
- Conventional Commits. Do not add dependencies on GWDistributions or any `GWBackground*`
  package, not even in the test target.
