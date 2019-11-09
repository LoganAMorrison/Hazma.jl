# include("theory.jl")
# include("constants.jl")
import Base.getproperty

abstract type AbstractScalarMediator <: AbstractTheoryMediator end

# mutable struct ScalarMediator <: AbstractScalarMediator
#     mχ::Real
#     ms::Real
#     gsχχ::Real
#     gsff::Real
#     gsGG::Real
#     gsFF::Real
#     Λ::Real
#
#     function ScalarMediator(mχ, ms, gsχχ, gsff, gsGG, gsFF, Λ)
#         if mχ < 0
#             throw(error("DM mass must be nonnegative"))
#         elseif ms < 0
#             throw(error("mediator mass must be nonnegative"))
#         elseif Λ < 0
#             throw(error("dimension 5 operator mass scale must be nonnegative"))
#         else
#             return new(mχ, ms, gsχχ, gsff, gsGG, gsFF, Λ)
#         end
#     end
# end

# mutable struct HiggsPortal <: AbstractScalarMediator
#     mχ::Real
#     ms::Real
#     gsχχ::Real
#     sinθ::Real
#
#     function HiggsPortal(mχ, ms, gsχχ, sinθ)
#         if mχ < 0
#             throw(error("DM mass must be nonnegative"))
#         elseif ms < 0
#             throw(error("mediator mass must be nonnegative"))
#         elseif sinθ < -1 || sinθ > 1
#             throw(error("sin(θ) must be between -1 and 1"))
#         else
#             return new(mχ, ms, gsχχ, sinθ)
#         end
#     end
# end

# mutable struct HeavyQuark <: AbstractScalarMediator
#     mχ::Real
#     ms::Real
#     gsχχ::Real
#     gsQQ::Real
#     mQ::Real
#     QQ::Real
#
#     function HeavyQuark(mχ, ms, gsχχ, gsQQ, mQ, QQ)
#         if mχ < 0
#             throw(error("DM mass must be nonnegative"))
#         elseif ms < 0
#             throw(error("mediator mass must be nonnegative"))
#         elseif mQ < 0
#             throw(error("heavy quark mass must be nonnegative"))
#         else
#             return new(mχ, ms, gsχχ, gsQQ, mQ, QQ)
#         end
#     end
# end

@make_mediator_theory ScalarMediator [
    (:mχ, nonnegative),
    (:ms, nonnegative),
    :gsχχ,
    :gsff,
    :gsGG,
    :gsFF,
    (:Λ, nonnegative),
] AbstractTheoryMediator
@make_mediator_theory HiggsPortal [
    (:mχ, nonnegative),
    (:ms, nonnegative),
    :gsχχ,
    (:sinθ, (sinθ) -> sinθ ≥ -1 && sinθ ≤ 1),
] AbstractScalarMediator
@make_mediator_theory HeavyQuark [
    (:mχ, nonnegative),
    (:ms, nonnegative),
    :gsχχ,
    :gsQQ,
    (:mQ, nonnegative),
    :QQ,
] AbstractScalarMediator

function Base.getproperty(mod::HiggsPortal, param::Symbol)
    if param == :gsff
        return mod.sinθ
    elseif param == :gsGG
        return 3 * mod.sinθ
    elseif param == :gsFF
        return -5 / 6 * mod.sinθ
    elseif param == :Λ
        return VH
    else
        return getfield(mod, param)
    end
end

function Base.getproperty(mod::HeavyQuark, param::Symbol)
    if param == :gsff
        return 0
    elseif param == :gsGG
        return mod.gsQQ
    elseif param == :gsFF
        return 2 * mod.gsQQ * mod.QQ^2
    elseif param == :Λ
        return mod.mQ
    else
        return getfield(mod, param)
    end
end

# Approximately true for all reasonable parameter values
scalar_vev(::AbstractScalarMediator) = 0.0
function list_decay_final_states(::AbstractScalarMediator)
    ["e⁺ e⁻", "μ⁺ μ⁻", "γ γ", "π⁰ π⁰", "π⁺ π⁻", "χ χ"]
end

function list_annihilation_final_states(::AbstractScalarMediator)
    ["e⁺ e⁻", "μ⁺ μ⁻", "γ γ", "π⁰ π⁰", "π⁺ π⁻", "s s"]
end

function list_decay_final_states(::HeavyQuark)
    ["γ γ", "π⁰ π⁰", "π⁺ π⁻", "χ χ"]
end

function list_annihilation_final_states(::HeavyQuark)
    ["γ γ", "π⁰ π⁰", "π⁺ π⁻", "s s"]
end

# --------------------- #
# Mediator decay mod.Γ_med #
# --------------------- #

function Γ_med(mod::AbstractScalarMediator, fs::String)
    if fs == "γ γ"
        return Γ_s_to_γγ(mod)
    elseif fs == "π⁰ π⁰"
        return Γ_s_to_π⁰π⁰(mod)
    elseif fs == "π⁺ π⁻"
        return Γ_s_to_ππ(mod)
    elseif fs == "χ χ"
        return Γ_s_to_χχ(mod)
    elseif fs == "e⁺ e⁻"
        return Γ_s_to_ee(mod)
    elseif fs == "μ⁺ μ⁻"
        return Γ_s_to_μμ(mod)
    elseif fs == "all"
        Γs = Dict(fs => Γ_med(mod, fs) for fs in list_decay_final_states(mod))
        Γs["total"] = sum(values(Γs))
        return Γs
    else
        return 0.0
    end
end

function Γ_s_to_γγ(mod::AbstractScalarMediator)
    (αem^2 * mod.gsFF^2 * model.ms^3) / (64 * mod.Λ^2 * π^3)
end

function Γ_s_to_π⁰π⁰(mod::AbstractScalarMediator)
    (mod.ms < 2mπ⁰) && return 0.0
    vs = scalar_vev(mod)

    ((sqrt(-4 * mπ⁰^2 + mod.ms^2) *
      (-162 * mod.gsGG * mod.Λ^3 * (-2 * mπ⁰^2 + mod.ms^2) * VH^2 +
       B0 * (md + mu) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
       (-3 * mod.Λ * VH + 3 * gsff * mod.Λ * vs + 2 * mod.gsGG * VH * vs) *
       (2 * mod.gsGG * VH * (9 * mod.Λ - 4 * mod.gsGG * vs) +
        9 * mod.gsff * mod.Λ * (3 * mod.Λ + 4 * mod.gsGG * vs)))^2) /
     (209952 * mod.Λ^6 * mod.ms^2 * π * VH^4 *
      (9 * mod.Λ + 4 * mod.gsGG * vs)^2))
end

function Γ_s_to_ππ(mod::AbstractScalarMediator)
    (mod.ms < 2mπ) && return 0.0
    vs = scalar_vev(mod)

    ((sqrt(-4 * mπ^2 + mod.ms^2) *
      (-162 * mod.gsGG * mod.Λ^3 * (-2 * mπ^2 + mod.ms^2) * VH^2 +
       B0 * (md + mu) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
       (-3 * mod.Λ * VH + 3 * mod.gsff * mod.Λ * vs + 2 * mod.gsGG * VH * vs) *
       (2 * mod.gsGG * VH * (9 * mod.Λ - 4 * mod.gsGG * vs) +
        9 * mod.gsff * mod.Λ * (3 * mod.Λ + 4 * mod.gsGG * vs)))^2) /
     (104976 * mod.Λ^6 * mod.ms^2 * π * VH^4 *
      (9 * mod.Λ + 4 * mod.gsGG * vs)^2))
end

function Γ_s_to_χχ(mod::AbstractScalarMediator)
    (mod.ms < 2 * mod.mχ) && return 0.0

    (mod.gsχχ^2 * (mod.ms^2 - 4 * mod.mχ^2)^1.5) / (8 * mod.ms^2 * π)
end

function Γ_s_to_ff(mod::AbstractScalarMediator, mf::Real)
    mod.ms < 2mf && return 0.0
    gsll = mod.gsff * (mf / VH)

    (gsll^2 * (-4 * mf^2 + mod.ms^2)^1.5) / (8 * mod.ms^2 * π)
end

Γ_s_to_ee(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, me)
Γ_s_to_μμ(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, mμ)
Γ_s_to_ee(mod::HeavyQuark) = 0.0;
Γ_s_to_μμ(mod::HeavyQuark) = 0.0;
# ----------------------------------- #
# DM self-annihilation cross sections #
# ----------------------------------- #

function σ_χχ(e_cm::Real, mod::AbstractScalarMediator, fs::String)
    if fs == "e⁺ e⁻"
        return σ_χχ_to_ee(e_cm, mod)
    elseif fs == "μ⁺ μ⁻"
        return σ_χχ_to_μμ(e_cm, mod)
    elseif fs == "γ γ"
        return σ_χχ_to_γγ(e_cm, mod)
    elseif fs == "π⁰ π⁰"
        return σ_χχ_to_π⁰π⁰(e_cm, mod)
    elseif fs == "π⁺ π⁻"
        return σ_χχ_to_ππ(e_cm, mod)
    elseif fs == "s s"
        return σ_χχ_to_ss(e_cm, mod)
    elseif fs == "all"
        σs = Dict(fs => σ_χχ(e_cm, mod, fs) for fs in list_annihilation_final_states(mod))
        σs["total"] = sum(values(σs))
        return σs
    else
        return zero(typeof(e_cm))
    end
end

function σ_χχ_to_ff(e_cm::Real, mod::AbstractScalarMediator, mf::Real)
    (e_cm < 2 * mf || e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))

    (mod.gsff^2 * mod.gsχχ^2 * mf^2 * (-4 * ml^2 + e_cm^2)^1.5 *
     sqrt(-4 * mod.mχ^2 + e_cm^2)) /
    (16 * π * e_cm^2 * VH^2 * ((mod.ms^2 - e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))
end

σ_χχ_to_ee(e_cm::Real, mod::AbstractScalarMediator) = σ_χχ_to_ff(e_cm, mod, me)
σ_χχ_to_μμ(e_cm::Real, mod::AbstractScalarMediator) = σ_χχ_to_ff(e_cm, mod, mμ)
σ_χχ_to_ee(e_cm::Real, mod::HeavyQuark) = zero(e_cm)
σ_χχ_to_μμ(e_cm::Real, mod::HeavyQuark) = zero(e_cm)

function σ_χχ_to_γγ(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))
    Γs = mod.Γ_med

    (αem^2 * mod.gsFF^2 * mod.gsχχ^2 * e_cm^3 * sqrt(-4 * mod.mχ^2 + e_cm^2)) /
    (128 * mod.Λ^2 * π^3 * ((mod.ms^2 - e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))
end

function σ_χχ_to_π⁰π⁰(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2 * mπ⁰ || e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))

    vs = scalar_vev(mod)

    (mod.gsχχ^2 * sqrt(-4 * mπ⁰^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2) *
     (162 * mod.gsGG * mod.Λ^3 * (2 * mπ⁰^2 - e_cm^2) * VH^2 +
      B0 * (md + mu) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
      (-3 * mod.Λ * VH + 3 * mod.gsff * mod.Λ * vs + 2 * mod.gsGG * VH * vs) *
      (2 * mod.gsGG * VH * (9 * mod.Λ - 4 * mod.gsGG * vs) +
       9 * mod.gsff * mod.Λ * (3 * mod.Λ + 4 * mod.gsGG * vs)))^2) /
    (419904 * mod.Λ^6 * π * e_cm^2 * VH^4 * (9 * mod.Λ + 4 * mod.gsGG * vs)^2 *
     ((mod.ms^2 - e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))
end

function σ_χχ_to_ππ(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2 * mπ || e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))

    (mod.gsχχ^2 * sqrt(-4 * mπ^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2) *
     (162 * mod.gsGG * mod.Λ^3 * (2 * mπ^2 - e_cm^2) * VH^2 +
      B0 * (md + mu) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
      (-3 * mod.Λ * VH + 3 * mod.gsff * mod.Λ * vs + 2 * mod.gsGG * VH * vs) *
      (2 * mod.gsGG * VH * (9 * mod.Λ - 4 * mod.gsGG * vs) +
       9 * mod.gsff * mod.Λ * (3 * mod.Λ + 4 * mod.gsGG * vs)))^2) /
    (209952 * mod.Λ^6 * π * e_cm^2 * VH^4 * (9 * mod.Λ + 4 * mod.gsGG * vs)^2 *
     ((mod.ms^2 - e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))
end

function σ_χχ_to_ss(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2 * mod.ms || e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))

    -(mod.gsχχ^4 *
      ((2 * sqrt(-4 * mod.ms^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2) *
        (3 * mod.ms^4 - 16 * mod.ms^2 * mod.mχ^2 +
         2 * mod.mχ^2 * (8 * mod.mχ^2 + e_cm^2))) /
       (mod.ms^4 - 4 * mod.ms^2 * mod.mχ^2 + mod.mχ^2 * e_cm^2) +
       ((6 * mod.ms^4 - 32 * mod.mχ^4 + 16 * mod.mχ^2 * e_cm^2 + e_cm^4 -
         4 * mod.ms^2 * (4 * mod.mχ^2 + e_cm^2)) *
        log((-2 * mod.ms^2 + e_cm^2 +
             sqrt(-4 * mod.ms^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2))^2 /
            (2 * mod.ms^2 - e_cm^2 +
             sqrt(-4 * mod.ms^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2))^2)) /
       (2 * mod.ms^2 - e_cm^2))) / (64 * π * e_cm^2 * (-4 * mod.mχ^2 + e_cm^2))
end

# ----------- #
# FSR Spectra #
# ----------- #

function dndeᵧ_χχ_to_ffγ(
    eᵧ::Real,
    e_cm::Real,
    mod::AbstractScalarMediator,
    mf::Real,
)
    e = eᵧ / e_cm
    rf = mf / e_cm
    s = e_cm^2 - 2 * e_cm * eᵧ

    (e_cm < 2 * mf || e_cm < 2 * mod.mχ || s < 4 * mf^2 || e_cm^2 < s) &&
    return zero(typeof(eᵧ))
    return (αem * (2 * (-1 + 4 * rf^2) *
             sqrt((-1 + 2 * e) * (-1 + 2 * e + 4 * rf^2)) +
             4 * (1 + 2 * (-1 + e) * e - 6 * rf^2 + 8 * e * rf^2 + 8 * rf^4) *
             atanh(sqrt(1 + (4 * rf^2) / (-1 + 2 * e))))) /
           (e * (1 - 4 * rf^2)^1.5 * π * e_cm)
end

dndeᵧ_χχ_to_eeγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator) =
    dndeᵧ_χχ_to_ffγ(eᵧ, e_cm, mod, me)
dndeᵧ_χχ_to_μμγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator) =
    dndeᵧ_χχ_to_ffγ(eᵧ, e_cm, mod, mμ)

# Zero these explicitly!
dndeᵧ_χχ_to_eeγ(eᵧ::Real, e_cm::Real, mod::HeavyQuark) = zero(eᵧ)
dndeᵧ_χχ_to_μμγ(eᵧ::Real, e_cm::Real, mod::HeavyQuark) = zero(eᵧ)

function dndeᵧ_χχ_to_ππγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)
    μπ = mπ / e_cm
    x = 2 * eᵧ / e_cm

    (x < 0 || 1 - 4 * μπ^2 < x || 2 * mod.mχ > e_cm) && return zero(typeof(eᵧ))
    dynamic = (-2 * sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x) +
               (-1 + 2 * μπ^2 + x) *
               log((1 - x - sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x))^2 /
                   (-1 + x - sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x))^2)) / x

    coeff = αem / (sqrt(1 - 4 * μπ^2) * π)

    2 * dynamic * coeff / e_cm
end

# ------------- #
# Decay spectra #
# ------------- #

dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator) =
    2 * dndeᵧ_μ_decay(eᵧ, e_cm / 2)

dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::HeavyQuark) = zero(eᵧ)

# TODO: other final states!

# ------------- #
# e⁺/e⁻ spectra #
# ------------- #

function dndeₑ_χχ(eₑ::Real, e_cm::Real, mod::AbstractScalarMediator, fs::String)
    if fs == "π⁺ π⁻"
        return dndeₑ_χχ_to_ππ(eₑ, e_cm, mod)
    elseif fs == "μ μ"
        return dndeₑ_χχ_to_μμ(eₑ, e_cm, mod)
    elseif fs == "s s"
        return dndeₑ_χχ_to_ss(eₑ, e_cm, mod)
    else
        return zero(typeof(eₑ))
    end
end

function lines_e(e_cm::T, mod::AbstractScalarMediator) where {T<:Real}
    Dict{String,T}("e⁺ e⁻" => e_cm / 2)
end

"""
    dndeₑ_ππ(eₑ::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into charged πons.
"""
dndeₑ_χχ_to_ππ(eₑ::Real, e_cm::Real, ::AbstractScalarMediator) =
    dndeₑ_π_decay(eₑ, e_cm / 2)

"""
    dndeₑ_μμ(eₑ::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into muons.
"""
dndeₑ_χχ_to_μμ(eₑ::Real, e_cm::Real, ::AbstractScalarMediator) =
    dndeₑ_μ_decay(eₑ, e_cm / 2)


"""
    dndeₑ_ss(eₑ::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into scalar mediators.
"""
function dndeₑ_χχ_to_ss(
    eₑ::Real,
    e_cm::Real,
    mod::AbstractScalarMediator;
    fs::String = "total",
)
    (e_cm < 2 * mod.mχ || e_cm < 2 * mod.ms) && return zero(typeof(eₑ))
    mod.Γ_med == 0 && return zero(typeof(eₑ))

    es = e_cm / 2

    # line contribution from s → e⁺ e⁻: factor of 2 for 2 electrons per S
    Γee = Γ_s_to_ee(mod) / mod.Γ_med
    lines_contrib = 2 * Γee *
                    boosted_delta_function(es, mod.ms, eₑ, me, mod.ms / 2)

    fs == "e⁺ e⁻" && return lines_contrib
    Γππ = Γ_s_to_ππ(mod) / mod.Γ_med
    Γμμ = Γ_s_to_μμ(mod) / mod.Γ_med
    # factors of 2 for two π (μ) per S
    dndeₑ_π_decay_srf(eₑ_srf) = 2 * Γππ * dndeₑ_π_decay(eₑ_srf, mod.ms / 2)
    dndeₑ_μ_decay_srf(eₑ_srf) = 2 * Γμμ * dndeₑ_μ_decay(eₑ_srf, mod.ms / 2)
    spectrum_rf(eₑ_srf) = dndeₑ_π_decay_srf(eₑ_srf) + dndeₑ_μ_decay_srf(eₑ_srf)

    if fs == "total"
        result = boost_spectrum(spectrum_rf, es, mod.ms, eₑ, me)
    elseif fs == "π⁺ π⁻"
        result = boost_spectrum(dndeₑ_π_decay_srf, es, mod.ms, eₑ, me)
    elseif fs == "μ⁺ μ⁻"
        result = boost_spectrum(dndeₑ_π_decay_srf, es, mod.ms, eₑ, me)
    else
        result = zero(typeof(eₑ))
    end

    # factor of 2 for 2 final state S's
    2 * (result + lines_contrib)
end
