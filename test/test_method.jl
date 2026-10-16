using BackgroundCosmology
using BackgroundCosmology: GaussLegendreRule
using ForwardDiff
using Integrals
using Test

# A third-party rule: composite trapezoid, defined only through `integrate` dispatch.
struct Trapezoid
    n::Int
end
function BackgroundCosmology.integrate(f, lo, hi, m::Trapezoid)
    h = (hi - lo) / m.n
    acc = (f(lo) + f(hi)) / 2
    for i in 1:(m.n - 1)
        acc += f(lo + i * h)
    end
    return h * acc
end

const COSMOLOGIES = (
    FlatLCDM(67.0, 0.315), FlatWCDM(67.0, 0.315, -0.9), FlatW0WaCDM(67.0, 0.315, -0.9, 0.2),
)
const _SCALAR_DEFAULT = BackgroundCosmology._SCALAR_RULE
const TIGHT = (; reltol = 1.0e-13, abstol = 0.0)

@testset "integrate" begin
    @test integrate(exp, 0.0, 1.0, GaussLegendreRule(10)) ≈ exp(1) - 1 rtol = 1.0e-14
    @test GaussLegendreRule(4).nodes isa Vector{Float64}
    # f's output type is kept, so dual-valued integrands work.
    @test integrate(x -> ForwardDiff.Dual(x, 1.0), 0.0, 1.0, GaussLegendreRule(3)) isa
        ForwardDiff.Dual
    err = try
        integrate(identity, 0.0, 1.0, :foo)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("Symbol", err.msg)
end

@testset "methods agree with the default" begin
    zs = (0.1, 1.0, 5.0, 20.0)
    for c in COSMOLOGIES, z in zs
        d0 = comoving_radial_distance(c, z)
        @test comoving_radial_distance(c, z; method = GaussLegendreRule(80)) ≈ d0 rtol = 1.0e-12
        @test comoving_radial_distance(c, z; method = (QuadGKJL(), TIGHT)) ≈ d0 rtol = 1.0e-11
        @test comoving_radial_distance(c, z; method = Integrals.GaussLegendre(; n = 60)) ≈ d0 rtol =
            1.0e-12
        # The composite Integrals rule goes through the same kernel.
        @test comoving_radial_distance(
            c, z; method = Integrals.GaussLegendre(; n = 20, subintervals = 8)
        ) ≈ d0 rtol = 1.0e-12
        @test comoving_radial_distance(c, z; method = Trapezoid(4000)) ≈ d0 rtol = 1.0e-6
    end
    for c in COSMOLOGIES, z in (0.5, 2.0, 1000.0)
        @test age(c, z; method = GaussLegendreRule(100)) ≈ age(c, z) rtol = 1.0e-12
        @test age(c, z; method = (QuadGKJL(), TIGHT)) ≈ age(c, z) rtol = 1.0e-11
    end
    c = COSMOLOGIES[3]
    z_grid = collect(LinRange(0.0, 20.0, 129))
    g0 = distance_and_volume_grid(c, z_grid)
    g = distance_and_volume_grid(c, z_grid; method = (QuadGKJL(), TIGHT))
    @test g.comoving_distance ≈ g0.comoving_distance rtol = 1.0e-11
    @test g.differential_comoving_volume ≈ g0.differential_comoving_volume rtol = 1.0e-11
end

@testset "wrappers forward `method`" begin
    c = COSMOLOGIES[2]
    m = Trapezoid(50) # deliberately crude, so the default can't pass by accident
    d = comoving_radial_distance(c, 2.0; method = m)
    @test d != comoving_radial_distance(c, 2.0)
    @test comoving_transverse_distance(c, 2.0; method = m) == d
    @test luminosity_distance(c, 2.0; method = m) == 3d
    @test angular_diameter_distance(c, 2.0; method = m) == d / 3
    @test comoving_volume(c, 2.0; method = m) == 4π / 3 * d^3
    @test distmod(c, 2.0; method = m) == 5 * log10(3d) + 25
    @test comoving_volume_element(c, 2.0; method = m) ==
        hubble_distance0(c) * d^2 / E(c, 2.0)
    # `age`'s integrand is a removable 0·∞ at √a = 0, so it needs an open rule.
    r = GaussLegendreRule(3)
    @test lookback_time(c, 1.5; method = r) == age(c, 0; method = r) - age(c, 1.5; method = r)
    @test lookback_time(c, 1.5; method = r) != lookback_time(c, 1.5)
end

@testset "ForwardDiff through other methods" begin
    build(θ) = FlatW0WaCDM(θ[1], θ[2], θ[3], θ[4])
    θ0 = [67.0, 0.315, -0.9, 0.2]
    quantities = (
        (θ, method) -> comoving_radial_distance(build(θ), 1.3; method),
        (θ, method) -> age(build(θ), 0.5; method),
    )
    methods = (GaussLegendreRule(80), (QuadGKJL(), TIGHT), Integrals.GaussLegendre(; n = 60))
    for q in quantities
        g0 = ForwardDiff.gradient(θ -> q(θ, _SCALAR_DEFAULT), θ0)
        for method in methods
            @test ForwardDiff.gradient(θ -> q(θ, method), θ0) ≈ g0 rtol = 1.0e-8
        end
    end
end

@testset "default path does not allocate" begin
    c = COSMOLOGIES[3]
    comoving_radial_distance(c, 1.3); age(c, 0.5)
    @test (@allocated comoving_radial_distance(c, 1.3)) == 0
    @test (@allocated age(c, 0.5)) == 0
end
