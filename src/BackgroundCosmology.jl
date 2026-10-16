module BackgroundCosmology

using FastGaussQuadrature: gausslegendre

export AbstractCosmology, FlatLCDM, FlatWCDM, FlatW0WaCDM,
    cosmology,
    H0, Ωm, E, H, w,
    hubble_dist0, hubble_dist,
    comoving_radial_dist, comoving_transverse_dist,
    angular_diameter_dist, luminosity_dist, distmod,
    comoving_volume, comoving_volume_element,
    age, lookback_time,
    hubble_constant_si,
    distance_and_volume_grid

include("conversion.jl")
include("model.jl")
include("quadrature.jl")
include("distance.jl")
include("time.jl")

end # module
