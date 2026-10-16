using BackgroundCosmology
using Reactant
using Test

Reactant.set_default_backend("cpu")

const MODELS = (
    (H0 = 67.0, Ωm = 0.315),
    (H0 = 67.0, Ωm = 0.315, w0 = -0.9),
    (H0 = 67.0, Ωm = 0.315, w0 = -0.9, wa = 0.2),
)
const QUANTITIES = (
    E, hubble_parameter, comoving_radial_distance, comoving_transverse_distance,
    angular_diameter_distance, luminosity_distance, distmod, comoving_volume,
    comoving_volume_element, age, lookback_time,
)

scalar_quantities(c, z) = map(f -> f(c, z), QUANTITIES)
changed_parameters(h) = map(x -> 0.95 * x, h)
host_grid(c) = distance_and_volume_grid(c, [0.0, 0.1, 0.2])
function constructor_values(C, h, z)
    extra = (; h..., ignored = 1.0)
    return (
        E(cosmology(; h...), z), E(cosmology(h), z),
        E(cosmology(C, extra), z), E(C(extra), z),
    )
end

function test_scalars(actual, expected)
    for (a, b) in zip(actual, expected)
        @test isapprox(Float64(a), b; rtol = 1.0e-12, atol = 0.0)
    end
    return
end

function test_grid(actual, expected)
    @test keys(actual) == keys(expected)
    for key in keys(expected)
        @test all(isapprox.(Array(actual[key]), expected[key]; rtol = 1.0e-12, atol = 0.0))
    end
    return
end

@testset "Reactant (Julia $VERSION, Reactant $(pkgversion(Reactant)))" begin
    @testset "constructors and scalar quantities" begin
        for h in MODELS
            c = cosmology(h)
            C = Base.typename(typeof(c)).wrapper
            h_ra = Reactant.to_rarray(h; track_numbers = true)
            c_ra = Reactant.to_rarray(c; track_numbers = true)
            @test cosmology(; h_ra...) isa C
            z_ra = Reactant.to_rarray(0.5; track_numbers = true)
            constructors = @compile constructor_values(C, h_ra, z_ra)
            test_scalars(constructors(C, h_ra, z_ra), constructor_values(C, h, 0.5))

            # Plain Julia redshifts are static; wrapped numbers vary at runtime.
            static_z = @compile scalar_quantities(c_ra, 0.5)
            test_scalars(static_z(c_ra, 0.5), scalar_quantities(c, 0.5))
            dynamic_z = @compile scalar_quantities(c_ra, z_ra)
            for z in (0.0, 0.1, 0.5, 2.0, 20.0)
                input_z = Reactant.to_rarray(z; track_numbers = true)
                test_scalars(dynamic_z(c_ra, input_z), scalar_quantities(c, z))
            end

            changed_h = changed_parameters(h)
            changed_c = cosmology(changed_h)
            changed_ra = Reactant.to_rarray(changed_c; track_numbers = true)
            changed_h_ra = Reactant.to_rarray(changed_h; track_numbers = true)
            test_scalars(dynamic_z(changed_ra, z_ra), scalar_quantities(changed_c, 0.5))
            test_scalars(static_z(changed_ra, 0.5), scalar_quantities(changed_c, 0.5))
            test_scalars(
                constructors(C, changed_h_ra, z_ra), constructor_values(C, changed_h, 0.5)
            )
        end
    end

    @testset "grid quantities" begin
        for h in MODELS
            c = cosmology(h)
            c_ra = Reactant.to_rarray(c; track_numbers = true)
            for z in ([0.0], [2.0], collect(range(0.0, 20.0; length = 1025)), [0.1, 0.2, 0.5])
                z_ra = Reactant.to_rarray(z)
                compiled = @compile distance_and_volume_grid(c_ra, z_ra)
                test_grid(compiled(c_ra, z_ra), distance_and_volume_grid(c, z))
                changed_z = z .* 0.9
                changed_c = cosmology(changed_parameters(h))
                changed_ra = Reactant.to_rarray(changed_c; track_numbers = true)
                test_grid(
                    compiled(changed_ra, Reactant.to_rarray(changed_z)),
                    distance_and_volume_grid(changed_c, changed_z)
                )
            end
            # A regular host grid remains accepted with traced cosmology parameters.
            compiled = @compile host_grid(c_ra)
            test_grid(compiled(c_ra), host_grid(c))
        end
    end

    @testset "grid tracing scales with array operations" begin
        c = Reactant.to_rarray(cosmology(); track_numbers = true)
        small = Reactant.to_rarray(collect(range(0.0, 2.0; length = 9)))
        large = Reactant.to_rarray(collect(range(0.0, 2.0; length = 1025)))
        small_ir = string(@code_hlo optimize = false distance_and_volume_grid(c, small))
        large_ir = string(@code_hlo optimize = false distance_and_volume_grid(c, large))
        small_ops = count("stablehlo.", small_ir)
        large_ops = count("stablehlo.", large_ir)
        @test small_ops > 0
        @test large_ops <= 2small_ops
    end
end
