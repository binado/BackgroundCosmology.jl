using ForwardDiff
using BackgroundCosmology
using Test

# Closed form for flat ΛCDM: t(a) = 2 / (3 H₀ √ΩΛ) · asinh(√(ΩΛ/Ωm) a^{3/2}).
function lcdm_age_gyr(H0, Ωm, z)
    ΩΛ = 1 - Ωm
    a = inv(1 + z)
    hubble_time_gyr = BackgroundCosmology._hubble_time_gyr(H0)
    return hubble_time_gyr * 2 / (3 * sqrt(ΩΛ)) * asinh(sqrt(ΩΛ / Ωm) * a^(3 / 2))
end

@testset "age matches the flat ΛCDM closed form" begin
    for (h0, om) in ((67.0, 0.315), (70.0, 0.3), (73.0, 0.25))
        c = FlatLCDM(h0, om)
        for z in (0.0, 0.1, 1.0, 5.0, 20.0, 1000.0)
            @test age(c, z) ≈ lcdm_age_gyr(h0, om, z) rtol = 1.0e-10
        end
    end
    # Sanity against the known number: ~13.8 Gyr for Planck-like parameters.
    @test 13.5 < age(cosmology()) < 14.0
    @test age(cosmology()) == age(cosmology(), 0.0)
end

@testset "lookback_time" begin
    for c in (
            FlatLCDM(67.0, 0.315), FlatWCDM(67.0, 0.315, -0.9),
            FlatW0WaCDM(67.0, 0.315, -0.9, 0.2),
        )
        @test lookback_time(c, 0.0) == 0.0
        @test lookback_time(c, 0) == 0.0
        zs = (0.1, 0.5, 1.0, 3.0)
        lb = [lookback_time(c, z) for z in zs]
        @test issorted(lb)
        for z in zs
            @test lookback_time(c, z) + age(c, z) ≈ age(c, 0) rtol = 1.0e-13
        end
        @test lookback_time(c, 20.0) < age(c, 0)
    end
end

@testset "age is differentiable" begin
    for (build, x0) in (
            (v -> FlatLCDM(67.0, v), 0.315),
            (v -> FlatLCDM(v, 0.315), 67.0),
            (v -> FlatWCDM(67.0, 0.315, v), -0.9),
            (v -> FlatW0WaCDM(67.0, 0.315, -0.9, v), 0.2),
        )
        f = v -> lookback_time(build(v), 1.5)
        d = ForwardDiff.derivative(f, x0)
        h = 1.0e-6
        @test isfinite(d)
        @test d ≈ (f(x0 + h) - f(x0 - h)) / (2h) rtol = 1.0e-6
    end
end
