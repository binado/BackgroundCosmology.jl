using Test
using BackgroundCosmology

@testset "BackgroundCosmology" begin
    include("test_cosmology.jl")
    include("test_time.jl")
    include("test_quality.jl")
end
