module BackgroundCosmologyIntegralsExt

using BackgroundCosmology: BackgroundCosmology, GaussLegendreRule
using Integrals: Integrals, IntegralProblem, IntegralFunction, solve
using Integrals.SciMLBase: AbstractIntegralAlgorithm, FullSpecialize

# `FullSpecialize` is essential: the default `AutoSpecialize` wraps the integrand in a
# type-erased callable, which costs ~4 µs and ~60 allocations per `solve` (about 7× the
# whole fixed-rule integral); specialized, the overhead is ~0.1 µs.
function BackgroundCosmology.integrate(
        f, lo, hi, alg::AbstractIntegralAlgorithm; kwargs...
    )
    g = IntegralFunction{false, FullSpecialize}((x, _) -> f(x), nothing)
    return solve(IntegralProblem(g, (lo, hi)), alg; kwargs...).u
end

# Tolerances, e.g. `(QuadGKJL(), (; reltol = 1.0e-12))`.
function BackgroundCosmology.integrate(
        f, lo, hi, (alg, kwargs)::Tuple{AbstractIntegralAlgorithm, NamedTuple}
    )
    return BackgroundCosmology.integrate(f, lo, hi, alg; kwargs...)
end

# Fixed Gauss–Legendre needs no problem/solve machinery: run our own kernel on its nodes,
# which is exactly as fast as the built-in `GaussLegendreRule`.
function BackgroundCosmology.integrate(f, lo, hi, alg::Integrals.GaussLegendre)
    rule = GaussLegendreRule(alg.nodes, alg.weights)
    h = (hi - lo) / alg.subintervals
    acc = BackgroundCosmology.integrate(f, lo, lo + h, rule)
    for i in 1:(alg.subintervals - 1)
        acc += BackgroundCosmology.integrate(f, lo + i * h, lo + (i + 1) * h, rule)
    end
    return acc
end

end # module
