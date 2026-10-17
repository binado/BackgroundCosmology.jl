using QuadGK
using Test
using ForwardDiff
using BackgroundCosmology
using BackgroundCosmology: hubble_constant_si, _de_density_ratio

@testset "hubble_constant_si" begin
    @test hubble_constant_si(70.0) ≈ 70.0 * 1000.0 / 3.085677581e22
end

@testset "basic cosmology helpers" begin
    c = FlatLCDM(67.0, 0.315)
    @test H0(c) == 67.0
    @test Ωm(c) == 0.315

    @test E(c, 0.0) ≈ 1.0
    @test hubble_parameter(c, 0.0) ≈ 67.0
    @test hubble_parameter(c, 1.0) ≈ 67.0 * E(c, 1.0)
    @test comoving_radial_distance(c, 0.0) ≈ 0.0

    # d_h = c / H0 in Mpc (speed of light in km/s over H0 in km/s/Mpc).
    @test hubble_distance(c) ≈ 299792.458 / 67.0

    z = [0.0, 0.1, 0.2]
    d_l = luminosity_distance.(c, z)
    @test d_l[1] ≈ 0.0
    @test d_l[3] > d_l[2] > d_l[1]
end

@testset "distance relations" begin
    for c in (
            FlatLCDM(67.0, 0.315), FlatWCDM(67.0, 0.315, -0.9),
            FlatW0WaCDM(67.0, 0.315, -0.9, 0.2),
        )
        for z in (0.1, 1.0, 3.0)
            d_c = comoving_radial_distance(c, z)
            @test comoving_transverse_distance(c, z) == d_c
            @test luminosity_distance(c, z) ≈ (1 + z) * d_c
            @test angular_diameter_distance(c, z) ≈ d_c / (1 + z)
            @test luminosity_distance(c, z) ≈ (1 + z)^2 * angular_diameter_distance(c, z)
            @test distmod(c, z) ≈ 5 * log10(luminosity_distance(c, z)) + 25
            @test comoving_volume(c, z) ≈ 4π / 3 * d_c^3
            @test comoving_volume_element(c, z) ≈ hubble_distance(c) * d_c^2 / E(c, z)
        end
    end

    # Hand-checkable low-redshift limit: d_c ≈ c z / H0 · (1 - (1 + q0) z / 2 ...).
    c = FlatLCDM(70.0, 0.3)
    @test luminosity_distance(c, 1.0e-4) ≈ 299792.458 / 70.0 * 1.0e-4 rtol = 1.0e-3

    # The volume element is the derivative of the volume out to z, per steradian.
    c = FlatLCDM(67.0, 0.315)
    z, h = 0.8, 1.0e-5
    dVdz = (comoving_volume(c, z + h) - comoving_volume(c, z - h)) / (2h)
    @test dVdz ≈ 4π * comoving_volume_element(c, z) rtol = 1.0e-7
end

@testset "cosmology construction" begin
    # Keyword constructor: the type follows the keywords that are given.
    @test cosmology(; H0 = 67.0, Ωm = 0.315) == FlatLCDM(67.0, 0.315)
    @test cosmology(; H0 = 67.0, Ωm = 0.315, w0 = -0.9) == FlatWCDM(67.0, 0.315, -0.9)
    @test cosmology(; H0 = 67.0, Ωm = 0.315, w0 = -0.9, wa = 0.2) ==
        FlatW0WaCDM(67.0, 0.315, -0.9, 0.2)
    @test_throws ArgumentError cosmology(; wa = 0.2)
    @test cosmology() isa FlatLCDM{Float64, Float64}
    @test @inferred(cosmology(; H0 = 67.0, Ωm = 0.315, w0 = -0.9)) isa FlatWCDM{Float64}

    # NamedTuple builders, used with sampler hyperparameters.
    h_lcdm = (H0 = 67.0, Ωm = 0.315)
    @test cosmology(FlatLCDM, h_lcdm) == FlatLCDM(67.0, 0.315)
    @test FlatLCDM(h_lcdm) == FlatLCDM(67.0, 0.315)
    @test cosmology(h_lcdm) == FlatLCDM(67.0, 0.315)
    # Extra propagation keys in `h` are ignored by the cosmology builder.
    @test cosmology(FlatLCDM, (; h_lcdm..., Ξ₀ = 1.2, Ξₙ = 2.0)) == FlatLCDM(67.0, 0.315)

    h_w0 = (; h_lcdm..., w0 = -0.9)
    @test cosmology(FlatWCDM, h_w0) == FlatWCDM(67.0, 0.315, -0.9)
    @test cosmology(h_w0) == FlatWCDM(67.0, 0.315, -0.9)

    h_cpl = (; h_w0..., wa = 0.2)
    @test cosmology(FlatW0WaCDM, h_cpl) == FlatW0WaCDM(67.0, 0.315, -0.9, 0.2)
    @test cosmology(h_cpl) == FlatW0WaCDM(67.0, 0.315, -0.9, 0.2)
end

@testset "dark-energy equation of state w" begin
    lcdm = FlatLCDM(67.0, 0.3)
    w0cdm = FlatWCDM(67.0, 0.3, -0.8)
    w0wacdm = FlatW0WaCDM(67.0, 0.3, -0.8, 0.3)
    for z in (0.0, 0.5, 1.0, 2.0)
        @test w(lcdm, z) ≈ -1.0
        @test w(w0cdm, z) ≈ -0.8
        @test w(w0wacdm, z) ≈ -0.8 + 0.3 * z / (1 + z)
    end
end

@testset "dark-energy density ratio" begin
    lcdm = FlatLCDM(67.0, 0.3)
    w0cdm = FlatWCDM(67.0, 0.3, -0.8)
    w0wacdm = FlatW0WaCDM(67.0, 0.3, -0.9, 0.2)
    for z in (0.1, 0.5, 1.0, 2.0, 5.0)
        @test _de_density_ratio(lcdm, z) ≈ 1.0
        # w0CDM closed form vs quadgk integral
        expected_w0, _ = quadgk(zp -> 3 * (1 + (-0.8)) / (1 + zp), 0.0, z; rtol = 1.0e-10)
        @test _de_density_ratio(w0cdm, z) ≈ exp(expected_w0) rtol = 1.0e-10
        # w0waCDM closed form vs quadgk integral
        expected_cpl, _ = quadgk(
            zp -> 3 * (1 + w(w0wacdm, zp)) / (1 + zp), 0.0, z; rtol = 1.0e-10
        )
        @test _de_density_ratio(w0wacdm, z) ≈ exp(expected_cpl) rtol = 1.0e-8
    end
end

@testset "comoving_radial_distance preserves AD tags at z=0" begin
    c_w0 = FlatWCDM(67.0, ForwardDiff.Dual(0.315), -0.9)
    @test comoving_radial_distance(c_w0, 0.0) ≈ 0.0
    @test comoving_radial_distance(c_w0, 0.0) isa ForwardDiff.Dual

    zs = [0.0, 0.1]
    r = comoving_radial_distance.(c_w0, zs)
    @test all(x -> x isa ForwardDiff.Dual, r)

    c_wa = FlatW0WaCDM(67.0, 0.315, -0.9, ForwardDiff.Dual(0.2))
    @test comoving_radial_distance(c_wa, 0.0) ≈ 0.0
    @test comoving_radial_distance(c_wa, 0.0) isa ForwardDiff.Dual
end

@testset "w preserves ForwardDiff derivatives" begin
    f_w0 = w0 -> w(FlatWCDM(67.0, 0.3, w0), 0.5)
    @test ForwardDiff.derivative(f_w0, -0.8) ≈ 1.0

    g_w0 = w0 -> E(FlatWCDM(67.0, 0.3, w0), 0.5)
    @test isfinite(ForwardDiff.derivative(g_w0, -0.8))
    @test ForwardDiff.derivative(g_w0, -0.8) != 0.0

    h_wa = wa -> w(FlatW0WaCDM(67.0, 0.3, -0.9, wa), 0.5)
    @test ForwardDiff.derivative(h_wa, 0.2) ≈ 0.5 / (1 + 0.5)
end

@testset "E(z) reduces to ΛCDM at w0=-1" begin
    lcdm = FlatLCDM(70.0, 0.3)
    w0cdm_lim = FlatWCDM(70.0, 0.3, -1.0)
    w0wacdm_lim = FlatW0WaCDM(70.0, 0.3, -1.0, 0.0)
    for z in (0.0, 0.1, 0.5, 2.0)
        @test E(lcdm, z) ≈ E(w0cdm_lim, z)
        @test E(lcdm, z) ≈ E(w0wacdm_lim, z)
    end
end

@testset "scalar distance path is pinned against adaptive QuadGK" begin
    # The fixed-order Gauss–Legendre rule replaced adaptive quadgk; order 40 is the
    # smallest order with margin below ~1e-12 relative up to z = 20 (32 misses it for
    # ΛCDM).
    for c in (
            FlatLCDM(67.0, 0.315), FlatWCDM(67.0, 0.315, -0.9),
            FlatW0WaCDM(67.0, 0.315, -0.9, 0.2),
        )
        for z in (0.1, 1.0, 2.5, 10.0, 20.0)
            reference, _ = quadgk(x -> inv(E(c, x)), 0.0, z; rtol = 1.0e-13)
            @test comoving_radial_distance(c, z) ≈ hubble_distance(c) * reference rtol = 1.0e-12
        end
    end
end

@testset "distance_and_volume_grid" begin
    z_grid = collect(LinRange(0.0, 20.0, 1025))

    # Grid and scalar paths are both Gauss–Legendre and agree to rounding error.
    for c in (
            FlatLCDM(67.0, 0.315), FlatWCDM(67.0, 0.315, -0.9),
            FlatW0WaCDM(67.0, 0.315, -0.9, 0.2),
        )
        g = distance_and_volume_grid(c, z_grid)
        @test g.comoving_distance ≈ comoving_radial_distance.(c, z_grid) rtol = 1.0e-12
        @test g.luminosity_distance ≈ luminosity_distance.(c, z_grid) rtol = 1.0e-12
        # The grid's volume is all-sky: 4π × the per-steradian element.
        @test g.differential_comoving_volume ≈
            4π .* comoving_volume_element.(c, z_grid) rtol = 1.0e-12
    end

    # astrogwb's contract: the grid need not start at zero — 0 is prepended
    # internally. A shifted grid agrees to rounding error, at the shared nodes, with the
    # zero-prepended grid: the leading zero-width interval contributes exactly 0
    # and every following interval is the same affine map. Prefix scans may
    # reassociate additions, so they need not be bit-identical.
    c = FlatLCDM(67.0, 0.315)
    z_shift = collect(LinRange(0.5, 20.0, 128))
    g_shift = distance_and_volume_grid(c, z_shift)
    g_prepended = distance_and_volume_grid(c, [0.0; z_shift])
    @test g_shift.comoving_distance ≈ g_prepended.comoving_distance[2:end] rtol = 1.0e-12
    @test g_shift.luminosity_distance ≈ g_prepended.luminosity_distance[2:end] rtol = 1.0e-12
    @test g_shift.differential_comoving_volume ≈
        g_prepended.differential_comoving_volume[2:end] rtol = 1.0e-12
    # A shifted grid at production density is also accurate, except for its first
    # interval [0, 0.5], which is ~25x wider than the rest and dominates the
    # composite-rule error (~1e-9 there). (The 128-point grid above is deliberately
    # coarse — validity is not accuracy, and the grid policy is caller-owned.)
    z_fine = collect(LinRange(0.5, 20.0, 1024))
    g_fine = distance_and_volume_grid(c, z_fine)
    @test g_fine.comoving_distance ≈ comoving_radial_distance.(c, z_fine) rtol = 1.0e-8

    # A single node is valid: one 4-node rule over [0, z]. Structurally sound, but
    # one wide interval is not an accurate composite rule — caller-owned policy.
    g_single = distance_and_volume_grid(c, [2.0])
    @test length(g_single.comoving_distance) == 1
    @test all(isfinite, g_single.comoving_distance)
    @test g_single.comoving_distance[1] > 0

    builders = (
        (v -> FlatLCDM(67.0, v), 0.315),
        (v -> FlatLCDM(v, 0.315), 67.0),
        (v -> FlatWCDM(67.0, 0.315, v), -0.9),
        (v -> FlatW0WaCDM(67.0, 0.315, -0.9, v), 0.2),
    )

    # ForwardDiff propagates through the full grid calculation.
    for (build, x0) in builders
        f_grid = v -> sum(distance_and_volume_grid(build(v), z_grid).luminosity_distance)
        d = ForwardDiff.derivative(f_grid, x0)
        @test isfinite(d)
        @test d != 0.0
        h = sqrt(eps(x0))
        @test d ≈ (f_grid(x0 + h) - f_grid(x0 - h)) / (2h) rtol = 1.0e-5
    end

    # And through the scalar path.
    for (build, x0) in builders
        f_scalar = v -> luminosity_distance(build(v), 1.5)
        d = ForwardDiff.derivative(f_scalar, x0)
        @test isfinite(d)
        @test d != 0.0
        h = sqrt(eps(x0))
        @test d ≈ (f_scalar(x0 + h) - f_scalar(x0 - h)) / (2h) rtol = 1.0e-5
    end

    g_zero = distance_and_volume_grid(c, [0.0])
    @test all(iszero, only.(values(g_zero)))
end

@testset "redshift derivatives at zero" begin
    for c in (FlatLCDM(67.0, 0.315), FlatWCDM(67.0, 0.315, -0.9), FlatW0WaCDM(67.0, 0.315, -0.9, 0.2))
        @test ForwardDiff.derivative(z -> comoving_radial_distance(c, z), 0.0) ≈ hubble_distance(c)
        @test ForwardDiff.derivative(z -> luminosity_distance(c, z), 0.0) ≈ hubble_distance(c)
    end
end
