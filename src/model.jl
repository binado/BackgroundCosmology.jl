"""
    AbstractCosmology

Abstract supertype for flat FLRW cosmologies without radiation. All quantities are plain
numbers: `H0` in km/s/Mpc, distances in Mpc, times in Gyr.
"""
abstract type AbstractCosmology end

"""
    FlatLCDM(H0, Ωm)

Flat ΛCDM cosmology (`w = -1`, radiation-free).
"""
struct FlatLCDM{TH0 <: Real, TΩm <: Real} <: AbstractCosmology
    H0::TH0
    Ωm::TΩm
end

"""
    FlatWCDM(H0, Ωm, w0)

Flat wCDM cosmology with constant dark-energy equation of state `w0`.
"""
struct FlatWCDM{TH0 <: Real, TΩm <: Real, Tw0 <: Real} <: AbstractCosmology
    H0::TH0
    Ωm::TΩm
    w0::Tw0
end

"""
    FlatW0WaCDM(H0, Ωm, w0, wa)

Flat w0waCDM (CPL) cosmology: `w(z) = w0 + wa·z/(1+z)`.
"""
struct FlatW0WaCDM{TH0 <: Real, TΩm <: Real, Tw0 <: Real, Twa <: Real} <: AbstractCosmology
    H0::TH0
    Ωm::TΩm
    w0::Tw0
    wa::Twa
end

"""    H0(c) -> Real

Present-day Hubble parameter in km/s/Mpc."""
H0(c::AbstractCosmology) = c.H0

"""    Ωm(c) -> Real

Present-day matter density parameter."""
Ωm(c::AbstractCosmology) = c.Ωm

Base.broadcastable(c::AbstractCosmology) = Ref(c)

"""
    cosmology(; H0 = 67.66, Ωm = 0.30966, w0 = nothing, wa = nothing)

Build a flat cosmology. The type follows which keywords are given: neither `w0` nor `wa`
→ [`FlatLCDM`](@ref); `w0` only → [`FlatWCDM`](@ref); both → [`FlatW0WaCDM`](@ref).
`wa` without `w0` is an error. The defaults are Planck 2018 (TT,TE,EE+lowE+lensing+BAO).

    cosmology(::Type{C}, h::NamedTuple)
    cosmology(h::NamedTuple)

Build from a hyperparameter `NamedTuple`; fields beyond those of `C` are ignored. Without
`C`, the type is inferred from the keys (`:wa` → `FlatW0WaCDM`, `:w0` → `FlatWCDM`, else
`FlatLCDM`).
"""
function cosmology(; H0 = 67.66, Ωm = 0.30966, w0 = nothing, wa = nothing)
    return _cosmology(H0, Ωm, w0, wa)
end

_cosmology(H0, Ωm, ::Nothing, ::Nothing) = FlatLCDM(H0, Ωm)
_cosmology(H0, Ωm, w0::Real, ::Nothing) = FlatWCDM(H0, Ωm, w0)
_cosmology(H0, Ωm, w0::Real, wa::Real) = FlatW0WaCDM(H0, Ωm, w0, wa)
function _cosmology(H0, Ωm, ::Nothing, ::Real)
    throw(ArgumentError("`wa` requires `w0` to be given"))
end

function cosmology(::Type{C}, h::NamedTuple) where {C <: AbstractCosmology}
    fn = fieldnames(C)
    return C(ntuple(i -> h[fn[i]], Val(length(fn)))...)
end

function cosmology(h::NamedTuple)
    :wa in keys(h) && return cosmology(FlatW0WaCDM, h)
    :w0 in keys(h) && return cosmology(FlatWCDM, h)
    return cosmology(FlatLCDM, h)
end

function (::Type{C})(h::NamedTuple) where {C <: AbstractCosmology}
    return cosmology(C, h)
end

"""
    w(c::AbstractCosmology, z) -> Real

Dark-energy equation of state `w(z)`.
"""
w(::FlatLCDM, z) = -one(z)
w(c::FlatWCDM, z) = c.w0
w(c::FlatW0WaCDM, z) = c.w0 + c.wa * z / (1 + z)

# ρ_DE(z) / ρ_DE(0): the closed-form integral of `w` through the Friedmann equation.
_de_density_ratio(::FlatLCDM, z) = one(z)
_de_density_ratio(c::FlatWCDM, z) = (1 + z)^(3 * (1 + c.w0))
function _de_density_ratio(c::FlatW0WaCDM, z)
    return (1 + z)^(3 * (1 + c.w0 + c.wa)) * exp(-3 * c.wa * z / (1 + z))
end

"""
    E(c::AbstractCosmology, z) -> Real

Dimensionless Hubble parameter `E(z) = H(z) / H₀`.
"""
function E(c::AbstractCosmology, z::Real)
    return sqrt(Ωm(c) * (1 + z)^3 + (1 - Ωm(c)) * _de_density_ratio(c, z))
end

"""
    H(c::AbstractCosmology, z) -> Real

Hubble parameter `H(z)` in km/s/Mpc.
"""
H(c::AbstractCosmology, z::Real) = H0(c) * E(c, z)
