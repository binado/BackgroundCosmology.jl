module BackgroundCosmology

using FastGaussQuadrature: gausslegendre

export AbstractCosmology, FlatLCDM, FlatWCDM, FlatW0WaCDM,
    cosmology,
    H0, Ωm, E, H, w,
    hubble_distance0, hubble_distance,
    comoving_radial_distance, comoving_transverse_distance,
    angular_diameter_distance, luminosity_distance, distmod,
    comoving_volume, comoving_volume_element,
    age, lookback_time,
    hubble_constant_si,
    distance_and_volume_grid,
    integrate

public GaussLegendreRule

include("conversion.jl")
include("model.jl")
include("quadrature.jl")
include("distance.jl")
include("time.jl")

end # module
