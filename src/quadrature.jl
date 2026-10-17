"""
    GaussLegendreRule(n::Integer)
    GaussLegendreRule(nodes, weights)

Fixed `n`-node Gauss–Legendre rule, with `nodes` and `weights` on the reference interval
`[-1, 1]`. This is the default `method` of every integral in the package; pass a different
order as `method = GaussLegendreRule(80)`.
"""
struct GaussLegendreRule{V <: AbstractVector{<:Number}}
    nodes::V
    weights::V
end

GaussLegendreRule(n::Integer) = GaussLegendreRule(gausslegendre(n)...)

"""
    integrate(f, lo, hi, method) -> typeof(f(midpoint))

Integral of the smooth scalar function `f` over `[lo, hi]` by `method`. Every integral in
the package goes through this function, with the change of variables and the grid
structure handled by the caller, so `method` only has to integrate a smooth one-dimensional
integrand.

Built-in methods are [`GaussLegendreRule`](@ref) and, with Integrals.jl loaded, any
SciML 1-D integral algorithm such as `QuadGKJL()` (give tolerances as
`(QuadGKJL(), (; reltol = 1e-12))`). To plug in your own rule, add a method for your type:

    BackgroundCosmology.integrate(f, lo, hi, ::MyRule) = ...

The result must follow `f`'s output type, so `ForwardDiff.Dual` cosmology parameters
propagate; fixed (non-adaptive) rules do this for free. The age integrand is a removable
`0·∞` at its lower endpoint, so `method` should be an open rule (Gauss–Legendre and
Gauss–Kronrod are; endpoint-evaluating rules such as trapezoid or Simpson return `NaN` for
[`age`](@ref)).
"""
function integrate(f, lo, hi, method)
    throw(
        ArgumentError(
            "no `integrate` method for $(typeof(method)); use a `GaussLegendreRule`, load " *
                "Integrals.jl to use its algorithms, or define " *
                "`BackgroundCosmology.integrate(f, lo, hi, ::$(typeof(method)))`"
        )
    )
end

function integrate(f, lo::Number, hi::Number, rule::GaussLegendreRule)
    half = (hi - lo) / 2
    mid = (lo + hi) / 2
    acc = zero(half * f(mid))
    @inbounds for k in eachindex(rule.nodes, rule.weights)
        acc += rule.weights[k] * f(mid + half * rule.nodes[k])
    end
    return half * acc
end

# Defaults, one per integrand (each is pinned against adaptive QuadGK in the tests).
# The grid path matches astrogwb's composite rule (`GAUSS_LEGENDRE_ORDER = 4`: 4 nodes
# per grid interval reach near machine precision for the smooth integrand 1/E(z)). The
# scalar rule covers [0, z] in one shot; order 40 is the smallest order with margin below
# 1e-12 relative error up to z = 20 for ΛCDM and w0waCDM. The age integral is smooth in
# √a, where order 64 is ample.
const _GRID_RULE = GaussLegendreRule(4)
const _SCALAR_RULE = GaussLegendreRule(40)
const _AGE_RULE = GaussLegendreRule(64)
