# Reference Gauss–Legendre nodes and weights on [-1, 1], computed once per rule.
# The grid path matches astrogwb's composite rule (`GAUSS_LEGENDRE_ORDER = 4`: 4 nodes
# per grid interval reach near machine precision for the smooth integrand 1/E(z)). The
# scalar rule covers [0, z] in one shot; order 40 is the smallest order with margin below
# 1e-12 relative error up to z = 20 for ΛCDM and w0waCDM (pinned against adaptive
# quadrature in the tests). The age integral is smooth in √a, where order 64 is ample.
const _GRID_GAUSS_LEGENDRE = gausslegendre(4)
const _SCALAR_GAUSS_LEGENDRE = gausslegendre(40)
const _AGE_GAUSS_LEGENDRE = gausslegendre(64)

"""
    _gauss_legendre_integral(f, a, b, (nodes, weights)) -> typeof(f(midpoint))

Integral of `f` over `[a, b]` with the reference Gauss–Legendre `nodes` and `weights` on
[-1, 1], mapped affinely onto the interval. The accumulation starts from `f`'s own output
type, so `ForwardDiff.Dual` cosmology parameters propagate.
"""
function _gauss_legendre_integral(f, a::Real, b::Real, (nodes, weights))
    half = (b - a) / 2
    mid = (a + b) / 2
    acc = zero(half * f(mid))
    @inbounds for k in eachindex(nodes, weights)
        acc += weights[k] * f(mid + half * nodes[k])
    end
    return half * acc
end
