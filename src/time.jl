# ∫₀^{√a} 2 ds / (s E(1/s² - 1)): the cosmic age in units of 1/H₀, with a = 1/(1+z) and
# s = √a. The integrand behaves like s²/√Ωm as s → 0 (matter domination), so unlike the
# 1/(a E) form in `a` it is smooth at the lower end and a fixed Gauss–Legendre rule
# converges geometrically.
function _age_integral(c::AbstractCosmology, z::Real, method)
    s_max = sqrt(inv(1 + z))
    g = s -> 2 * inv(s * E(c, inv(s^2) - 1))
    return integrate(g, zero(s_max), s_max, method)
end

"""
    age(c, z = 0) -> Real

Age of the universe at redshift `z` in Gyr, `(1/H₀) ∫₀ᵃ da′ / (a′ E(a′))` with
`a = 1/(1+z)`. Integrated in `√a`, which removes the matter-era endpoint singularity
(radiation is not modeled), by `method` (see [`integrate`](@ref)); the default is a fixed
64-node Gauss–Legendre rule.
"""
function age(c::AbstractCosmology, z::Real = 0; method = _AGE_RULE)
    return _hubble_time_gyr(H0(c)) * _age_integral(c, z, method)
end

"""
    lookback_time(c, z) -> Real

Time elapsed between redshift `z` and today, `age(c, 0) - age(c, z)`, in Gyr.
"""
function lookback_time(c::AbstractCosmology, z::Real; kwargs...)
    return age(c, 0; kwargs...) - age(c, z; kwargs...)
end
