"""Speed of light in km/s."""
const SPEED_OF_LIGHT_KM_S = 299792.458

"""1 Mpc in meters (IAU/CODATA convention)."""
const METERS_PER_MPC = 3.085677581e22

"""1 Gyr in seconds (Julian years)."""
const SECONDS_PER_GYR = 365.25 * 24 * 3600 * 1.0e9

"""
    hubble_constant_si(H0_km_s_mpc::Real) -> Float64

Hubble constant in s⁻¹ from ``H_0`` in **km/s/Mpc** (same units as the `H0` of every
[`AbstractCosmology`](@ref)).
"""
function hubble_constant_si(H0_km_s_mpc::Real)
    return Float64(H0_km_s_mpc) * 1000.0 / METERS_PER_MPC
end

# 1 / H0 in Gyr for H0 in km/s/Mpc.
# Deliberately no `Float64` conversion, so ForwardDiff duals in `H0` propagate.
_hubble_time_gyr(H0_km_s_mpc::Real) = METERS_PER_MPC / (1000 * H0_km_s_mpc) / SECONDS_PER_GYR
