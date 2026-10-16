# ∫₀^{√a} 2 ds / (s E(1/s² - 1)): the cosmic age in units of 1/H₀, with a = 1/(1+z) and
# s = √a. The integrand behaves like s²/√Ωm as s → 0 (matter domination), so unlike the
# 1/(a E) form in `a` it is smooth at the lower end and a fixed Gauss–Legendre rule
# converges geometrically.
function _age_integral(c::AbstractCosmology, z::Real)
    s_max = sqrt(inv(1 + z))
    g = s -> 2 * inv(s * E(c, inv(s^2) - 1))
    return _gauss_legendre_integral(g, zero(s_max), s_max, _AGE_GAUSS_LEGENDRE)
end

"""
    age(c, z = 0) -> Real

Age of the universe at redshift `z` in Gyr, `(1/H₀) ∫₀ᵃ da′ / (a′ E(a′))` with
`a = 1/(1+z)`. Integrated by a fixed 64-node Gauss–Legendre rule in `√a`, which removes
the matter-era endpoint singularity (radiation is not modeled).
"""
age(c::AbstractCosmology, z::Real = 0) = _hubble_time_gyr(H0(c)) * _age_integral(c, z)

"""
    lookback_time(c, z) -> Real

Time elapsed between redshift `z` and today, `age(c, 0) - age(c, z)`, in Gyr.
"""
lookback_time(c::AbstractCosmology, z::Real) = age(c, 0) - age(c, z)
