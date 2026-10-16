"""
    hubble_distance(c) -> Real

Hubble distance `c / H₀` in Mpc.
"""
hubble_distance(c::AbstractCosmology) = SPEED_OF_LIGHT_KM_S / H0(c)

"""
    hubble_distance(c, z) -> Real

Hubble distance at redshift `z`, `c / H(z)`, in Mpc.
"""
hubble_distance(c::AbstractCosmology, z::Real) = hubble_distance(c) / E(c, z)

"""
    comoving_radial_distance(c, z) -> Real

Line-of-sight comoving distance `d_c(z) = (c/H₀) ∫₀ᶻ dz′/E(z′)` in Mpc. The integral is
evaluated by `method` (see [`integrate`](@ref)) over `[0, z]`; the default is one fixed
40-node Gauss–Legendre rule (no adaptive branching, so ForwardDiff duals propagate).

All the distance and volume functions below forward `method` to this integral.
"""
function comoving_radial_distance(
        c::AbstractCosmology, z::Real; method = _SCALAR_RULE
    )
    d_h = hubble_distance(c)
    iszero(z) && return zero(d_h * z * inv(E(c, z)))
    integral = integrate(x -> inv(E(c, x)), zero(z), z, method)
    return d_h * integral
end

"""
    comoving_transverse_distance(c, z) -> Real

Transverse comoving distance in Mpc. The models are flat, so this is
[`comoving_radial_distance`](@ref).
"""
function comoving_transverse_distance(c::AbstractCosmology, z::Real; kwargs...)
    return comoving_radial_distance(c, z; kwargs...)
end

"""
    angular_diameter_distance(c, z) -> Real

Angular-diameter distance `d_c / (1 + z)` in Mpc.
"""
function angular_diameter_distance(c::AbstractCosmology, z::Real; kwargs...)
    return comoving_transverse_distance(c, z; kwargs...) / (1 + z)
end

"""
    luminosity_distance(c, z) -> Real

Luminosity distance `(1 + z) d_c` in Mpc.
"""
function luminosity_distance(c::AbstractCosmology, z::Real; kwargs...)
    return (1 + z) * comoving_transverse_distance(c, z; kwargs...)
end

"""
    distmod(c, z) -> Real

Distance modulus `5 log₁₀(d_L / Mpc) + 25`.
"""
function distmod(c::AbstractCosmology, z::Real; kwargs...)
    return 5 * log10(luminosity_distance(c, z; kwargs...)) + 25
end

"""
    comoving_volume(c, z) -> Real

Full-sky comoving volume `4π d_c³ / 3` out to redshift `z`, in Mpc³.
"""
function comoving_volume(c::AbstractCosmology, z::Real; kwargs...)
    return 4π / 3 * comoving_transverse_distance(c, z; kwargs...)^3
end

"""
    comoving_volume_element(c, z) -> Real

Comoving volume element **per steradian**, `dV_c/dz/dΩ = (c/H₀) d_c² / E(z)`, in
Mpc³ sr⁻¹ (same convention as Cosmology.jl). Multiply by `4π` for the all-sky value, which
is what [`distance_and_volume_grid`](@ref) tabulates.
"""
function comoving_volume_element(c::AbstractCosmology, z::Real; kwargs...)
    return hubble_distance(c) * comoving_transverse_distance(c, z; kwargs...)^2 / E(c, z)
end

"""
    distance_and_volume_grid(c::AbstractCosmology, z; method = GaussLegendreRule(4))
        -> (; comoving_distance, luminosity_distance, differential_comoving_volume)

Tabulate three distance quantities on `z` in a single pass, sharing one `1/E(z)`
evaluation pass and one cumulative quadrature between them. Distances are in Mpc.

`differential_comoving_volume` is the **all-sky** `dV_c/dz = 4π × `
[`comoving_volume_element`](@ref), in Mpc³ per unit redshift — the `4π` lives here, not in
the redshift-distribution consumers, matching the Python `astrogwb` stack.

The redshift integral is composite: `method` (default a fixed 4-node Gauss–Legendre
rule) is applied on each grid interval, with `0` prepended internally — matching astrogwb's `distance_and_volume_grid`
to rounding error for the smooth integrand `1/E(z)`: `d_c = d_h · cumsum(∫ 1/E)`,
`d_L = (1+z) · d_c`, `dV_c/dz = 4π · d_h · d_c² / E(z)`. This is the efficient batched
path for models that already evaluate and normalize quantities on a redshift grid; the
scalar functions use a single 40-node rule over `[0, z]` instead. The grid approximation
and any interpolation policy belong to the caller.

Takes the grid **array**, not `(z_min, z_max, n)`, so the caller's grid is the grid used
— there is no way for tabulation and interpolation to disagree about nodes. The grid must
be non-negative and strictly increasing, with at least one node (astrogwb's contract): a
grid that does not start at zero is fine, because the cumulative comoving-distance
integral always starts from the internally prepended `d_c(0) = 0`.
"""
function distance_and_volume_grid(
        c::AbstractCosmology, z::AbstractVector{<:Real}; method = _GRID_RULE
    )
    length(z) >= 1 || throw(
        ArgumentError(
            "distance_and_volume_grid requires at least one grid point"
        )
    )
    all(>=(0), z) || throw(
        ArgumentError(
            "distance_and_volume_grid requires a non-negative grid"
        )
    )
    all(>(0), diff(z)) || throw(
        ArgumentError(
            "distance_and_volume_grid requires a strictly increasing grid"
        )
    )
    inv_E = inv.(E.(c, z))
    # Cumulative ∫₀^z 1/E over each interval [z_{i-1}, z_i] (z_{-1} = 0 prepended),
    # accumulated left to right so the summation order matches astrogwb's `cumsum`.
    integral = similar(inv_E)
    acc = zero(eltype(integral))
    z_prev = zero(eltype(z))
    @inbounds for (i, z_i) in enumerate(z)
        acc += integrate(x -> inv(E(c, x)), z_prev, z_i, method)
        integral[i] = acc
        z_prev = z_i
    end
    d_h = hubble_distance(c)
    d_c = d_h .* integral
    return (;
        comoving_distance = d_c,
        luminosity_distance = (1 .+ z) .* d_c,
        differential_comoving_volume = @. 4π * d_h * d_c^2 * inv_E
    )
end
