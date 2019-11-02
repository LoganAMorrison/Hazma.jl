# include("theory.jl")
# include("constants.jl")
import Base.getproperty

abstract type AbstractScalarMediator <: AbstractTheoryMediator end

# Generic model
mutable struct ScalarMediator <: AbstractScalarMediator
    mχ::Real
    ms::Real
    gsχχ::Real
    gsff::Real
    gsGG::Real
    gsFF::Real
    Λ::Real

    function ScalarMediator(mχ, ms, gsχχ, gsff, gsGG, gsFF, Λ)
        if mχ < 0
            throw(error("DM mass must be nonnegative"))
        elseif ms < 0
            throw(error("mediator mass must be nonnegative"))
        elseif Λ < 0
            throw(error("dimension 5 operator mass scale must be nonnegative"))
        else
            return new(mχ, ms, gsχχ, gsff, gsGG, gsFF, Λ)
        end
    end
end

mutable struct HiggsPortal <: AbstractScalarMediator
    mχ::Real
    ms::Real
    gsχχ::Real
    sinθ::Real

    function HiggsPortal(mχ, ms, gsχχ, sinθ)
        if mχ < 0
            throw(error("DM mass must be nonnegative"))
        elseif ms < 0
            throw(error("mediator mass must be nonnegative"))
        elseif sinθ < -1 || sinθ > 1
            throw(error("sin(θ) must be between -1 and 1"))
        else
            return new(mχ, ms, gsχχ, sinθ)
        end
    end
end

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

mutable struct HeavyQuark <: AbstractScalarMediator
    mχ::Real
    ms::Real
    gsχχ::Real
    gsQQ::Real
    mQ::Real
    QQ::Real

    function HeavyQuark(mχ, ms, gsχχ, gsQQ, mQ, QQ)
        if mχ < 0
            throw(error("DM mass must be nonnegative"))
        elseif ms < 0
            throw(error("mediator mass must be nonnegative"))
        elseif mQ < 0
            throw(error("heavy quark mass must be nonnegative"))
        else
            return new(mχ, ms, gsχχ, gsQQ, mQ, QQ)
        end
    end
end

function Base.getproperty(mod::HeavyQuark, param::Symbol)
    if param == :gsff
        return 0.0
    elseif param == :gsGG
        return mod.gsQQ
    elseif param == :gsFF
        return 2.0 * mod.gsQQ * mod.QQ^2
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
# Mediator decay widths #
# --------------------- #

# TODO: generate with macro
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
    (αem^2 * mod.gsFF^2 * mod.ms^3) / (128.0 * mod.Λ^2 * π^3)
end

function Γ_s_to_π⁰π⁰(mod::AbstractScalarMediator)
    (mod.ms < 2mπ⁰) && return 0.0
    vs = scalar_vev(mod)

    return (sqrt(-4 * mπ⁰^2 + mod.ms^2) *
        ((mod.gsGG * (4 * mπ⁰^2 - 2 * mod.ms^2)) /
         (9 * mod.Λ + 4 * mod.gsGG * vs) +
         (B0 * (md + mu) *
          (27 * mod.gsff^2 * mod.Λ^2 * vs * (3 * mod.Λ + 4 * mod.gsGG * vs) -
           2 * mod.gsGG * VH^2 *
           (27 * mod.Λ^2 - 30 * mod.gsGG * mod.Λ * vs + 8 * mod.gsGG^2 * vs^2) +
           mod.gsff *
           (-81 * mod.Λ^3 * VH + 48 * mod.gsGG^2 * mod.Λ * VH * vs^2))) /
         (81.0 * mod.Λ^3 * VH^2))^2) / (16.0 * mod.ms^2 * π)
end

function Γ_s_to_ππ(mod::AbstractScalarMediator)
    (mod.ms < 2mπ) && return 0.0
    vs = scalar_vev(mod)
    return (sqrt(-4 * mπ^2 + mod.ms^2) *
            ((mod.gsGG * (4 * mπ^2 - 2 * mod.ms^2)) /
             (9 * mod.Λ + 4 * mod.gsGG * vs) +
             (B0 * (md + mu) *
              (27 * mod.gsff^2 * mod.Λ^2 * vs *
               (3 * mod.Λ + 4 * mod.gsGG * vs) - 2 * mod.gsGG * VH^2 *
               (27 * mod.Λ^2 - 30 * mod.gsGG * mod.Λ * vs +
                8 * mod.gsGG^2 * vs^2) +
               mod.gsff *
               (-81 * mod.Λ^3 * VH + 48 * mod.gsGG^2 * mod.Λ * VH * vs^2))) /
             (81.0 * mod.Λ^3 * VH^2))^2) / (16.0 * mod.ms^2 * pi)
end

function Γ_s_to_χχ(mod::AbstractScalarMediator)
    (mod.ms < 2 * mod.mχ) && return 0.0
    return (mod.gsχχ^2 * (mod.ms - 2 * mod.mχ) * (mod.ms + 2 * mod.mχ) *
            sqrt(mod.ms^2 - 4 * mod.mχ^2)) / (32.0 * mod.ms^2 * pi)
end

function Γ_s_to_ff(mod::AbstractScalarMediator, mf::Real)
    mod.ms < 2mf && return 0.0
    return (-(mod.gsff^2 * mf^2 * (2 * mf - mod.ms) * (2 * mf + mod.ms) *
              sqrt(-4 * mf^2 + mod.ms^2)) / (32.0 * mod.ms^2 * π * VH^2))
end

# TODO: generate with macro
Γ_s_to_ee(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, me)
Γ_s_to_μμ(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, mμ)
Γ_s_to_ee(mod::HeavyQuark) = 0.0
Γ_s_to_μμ(mod::HeavyQuark) = 0.0

# ----------------------------------- #
# DM self-annihilation cross sections #
# ----------------------------------- #

# TODO: generate with macro
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
        return 0.
    end
end

# TODO: WRONG!
function σ_χχ_to_ff(e_cm::Real, mod::AbstractScalarMediator, mf::Real)
    (e_cm < 2.0mf || e_cm < 2.0 * mod.mχ) && return 0.0
    Γs = Γ_med(mod)
    return ((mod.gsFF^2 * mod.gsχχ^2 * mf^2 * (-2.0 * mod.mχ + e_cm) *
             (2.0 * mod.mχ + e_cm) * (-4.0 * mf^2 + e_cm^2)^1.5) /
            (16.0π * e_cm^2 * sqrt(e_cm^2 - 4.0 * mod.mχ^2) * VH^2 *
             (mod.ms^4 - 2.0 * mod.ms^2 * e_cm^2 + e_cm^4 + mod.ms^2 * Γs^2)))
end

# TODO: generate with macro
σ_χχ_to_ee(e_cm::Real, mod::AbstractScalarMediator) =
    σ_χχ_to_ff(e_cm, mod, me)
σ_χχ_to_μμ(e_cm::Real, mod::AbstractScalarMediator) =
    σ_χχ_to_ff(e_cm, mod, mμ)
σ_χχ_to_ee(e_cm::Real, mod::HeavyQuark) = zero(e_cm)
σ_χχ_to_μμ(e_cm::Real, mod::HeavyQuark) = zero(e_cm)

function σ_χχ_to_γγ(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2.0 * mod.mχ) && return 0.0

    Γs = Γ_med(mod)
    rχs = mod.mχ / e_cm
    return ((αem^2 * mod.gsFF^2 * mod.gsχχ^2 * e_cm^4 * sqrt(1.0 - 4.0 * rχs^2)) /
            (64.0 * mod.Λ^2 * π^3 *
             (mod.ms^4 + e_cm^4 + mod.ms^2 * (-2.0 * e_cm^2 + Γs^2))))
end

function σ_χχ_to_π⁰π⁰(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2.0mπ⁰ || e_cm < 2.0 * mod.mχ) && return 0.0

    Γs = Γ_med(mod)
    vs = scalar_vev(mod)
    rπ⁰s = mπ⁰ / e_cm
    rχs = mod.mχ / e_cm

    return ((mod.gsχχ^2 * sqrt((-1.0 + 4.0 * rπ⁰s^2) * (-1.0 + 4.0 * rχs^2)) *
             (162.0 * mod.gsGG * mod.Λ^3 * e_cm^2 * (-1.0 + 2.0 * rπ⁰s^2) * VH^2 +
              B0 * (md + mu) * (9.0 * mod.Λ + 4.0 * mod.gsGG * vs) *
              (27.0 * mod.gsff^2 * mod.Λ^2 * vs *
               (3.0 * mod.Λ + 4.0 * mod.gsGG * vs) - 2.0 * mod.gsGG * VH^2 *
               (27.0 * mod.Λ^2 - 30.0 * mod.gsGG * mod.Λ * vs +
                8.0 * mod.gsGG^2 * vs^2) +
               mod.gsff *
               (-81.0 * mod.Λ^3 * VH + 48.0 * mod.gsGG^2 * mod.Λ * VH * vs^2)))^2) /
            (209952.0 * mod.Λ^6 * π * VH^4 * (9.0 * mod.Λ + 4.0 * mod.gsGG * vs)^2 *
             (mod.ms^4 + e_cm^4 + mod.ms^2 * (-2.0 * e_cm^2 + Γs^2))))
end

function σ_χχ_to_ππ(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2.0 * mπ || e_cm < 2.0 * mod.mχ) && return 0.0

    Γs = Γ_med(mod)
    vs = scalar_vev(mod)
    rπs = mπ / e_cm
    rχs = mod.mχ / e_cm

    return ((mod.gsχχ^2 * sqrt((-1.0 + 4.0 * rπs^2) * (-1.0 + 4.0 * rχs^2)) *
             (162.0 * mod.gsGG * mod.Λ^3 * e_cm^2 * (-1.0 + 2.0 * rπs^2) * VH^2 +
              B0 * (md + mu) * (9.0 * mod.Λ + 4.0 * mod.gsGG * vs) *
              (27.0 * mod.gsff^2 * mod.Λ^2 * vs *
               (3.0 * mod.Λ + 4.0 * mod.gsGG * vs) - 2.0 * mod.gsGG * VH^2 *
               (27.0 * mod.Λ^2 - 30.0 * mod.gsGG * mod.Λ * vs +
                8.0 * mod.gsGG^2 * vs^2) +
               mod.gsff *
               (-81.0 * mod.Λ^3 * VH + 48.0 * mod.gsGG^2 * mod.Λ * VH * vs^2)))^2) /
            (104976.0 * mod.Λ^6 * π * VH^4 * (9.0 * mod.Λ + 4.0 * mod.gsGG * vs)^2 *
             (mod.ms^4 + e_cm^4 + mod.ms^2 * (-2.0 * e_cm^2 + Γs^2))))
end

function σ_χχ_to_ss(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2.0 * mod.ms || e_cm < 2.0 * mod.mχ) && return 0.0

    return -((mod.gsχχ^4 * sqrt(-4.0 * mod.ms^2 + e_cm^2) *
              sqrt(-4.0 * mod.mχ^2 + e_cm^2) *
              (-2.0 / (4.0 * mod.mχ^2 - e_cm^2) -
               (mod.ms^2 - 4.0 * mod.mχ^2)^2 / ((4.0 * mod.mχ^2 - e_cm^2) *
                (mod.ms^4 - 4.0 * mod.ms^2 * mod.mχ^2 + mod.mχ^2 * e_cm^2)) -
               (2.0 *
                (6.0 * mod.ms^4 - 32.0 * mod.mχ^4 + 16.0 * mod.mχ^2 * e_cm^2 +
                 e_cm^4 - 4.0 * mod.ms^2 * (4.0 * mod.mχ^2 + e_cm^2)) *
                atanh((sqrt(-4.0 * mod.ms^2 + e_cm^2) * sqrt(-4.0 * mod.mχ^2 + e_cm^2)) /
                      (-2.0 * mod.ms^2 + e_cm^2))) /
               (sqrt(-4.0 * mod.ms^2 + e_cm^2) * (-2.0 * mod.ms^2 + e_cm^2) *
                (-4.0 * mod.mχ^2 + e_cm^2)^1.5))) / (16.0 * π * e_cm^2))
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
    s = e_cm^2 - 2.0 * e_cm * eᵧ

    (e_cm < 2.0mf || e_cm < 2.0 * mod.mχ || s < 4.0 * mf^2 || e_cm^2 < s) &&
    return 0.0

    return (αem * (2.0 * (-1.0 + 4.0 * rf^2) *
             sqrt((-1.0 + 2.0 * e) * (-1.0 + 2.0 * e + 4.0 * rf^2)) + 4.0 *
             (1.0 + 2.0 * (-1.0 + e) * e - 6.0 * rf^2 + 8.0 * e * rf^2 +
              8.0 * rf^4) *
             atanh(sqrt(1 + (4.0 * rf^2) / (-1.0 + 2.0 * e))))) /
           (e * (1.0 - 4.0 * rf^2)^1.5 * π * e_cm)
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
    x = 2.0eᵧ / e_cm

    (x < 0 || 1.0 - 4.0 * μπ^2 < x || 2.0 * mod.mχ > e_cm) && return 0.0

    dynamic = (-2.0 * sqrt(1.0 - x) * sqrt(1.0 - 4.0 * μπ^2 - x) +
               (-1.0 + 2.0 * μπ^2 + x) *
               log((1.0 - x - sqrt(1.0 - x) * sqrt(1.0 - 4.0 * μπ^2 - x))^2 /
                   (-1.0 + x - sqrt(1.0 - x) * sqrt(1.0 - 4.0 * μπ^2 - x))^2)) /
              x

    coeff = αem / (sqrt(1.0 - 4.0 * μπ^2) * π)

    2.0 * dynamic * coeff / e_cm
end

# ------------- #
# Decay spectra #
# ------------- #

function dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)
    return 2.0 * decay_spectrum_muon(eᵧ, e_cm / 2.0)
end

# TODO: generate with macro
dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::HeavyQuark) = zero(eᵧ)

# TODO: other final states!

# ------------- #
# e⁺/e⁻ spectra #
# ------------- #

# TODO: generate with a macro
function dndeₑ(eₑ::Real, e_cm::Real, mod::AbstractScalarMediator, fs::String)
    if fs == "π⁺ π⁻"
        return dndeₑ_ππ(eₑ, e_cm, mod)
    elseif fs == "μ μ"
        return dndeₑ_μμ(eₑ, e_cm, mod)
    elseif fs == "s s"
        return dndeₑ_ss(eₑ, e_cm, mod)
    else
        return 0.0
    end
end

function lines_e(e_cm::T, mod::AbstractScalarMediator) where T <: Real
    return Dict{String,T}(
        "e⁺ e⁻" => e_cm / 2.0
    )
end

"""
    dndeₑ_ππ(eₑ::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into charged pions.
"""
dndeₑ_ππ(eₑ::Real, e_cm::Real, ::AbstractScalarMediator) =
    positron_spectrum_charged_pion(eₑ, e_cm / 2.0)

"""
    dndeₑ_μμ(eₑ::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into muons.
"""
dndeₑ_μμ(eₑ::Real, e_cm::Real, ::AbstractScalarMediator) =
    positron_spectrum_muon(eₑ, e_cm / 2.0)

function dndeₑ_ss_integrand(
    cosθ::Real,
    eₑ::Real,
    es::Real,
    mod::AbstractScalarMediator,
    fs::String,
)
    eₑ < me && return 0.0

    pwμμ = Γ_s_to_ff(mod, mμ)
    pwππ = Γ_s_to_ππ(mod)

    p = sqrt(eₑ^2 - me^2)
    γ = es / ms
    β = sqrt(1.0 - (ms / es)^2)
    eₑ_srf = γ * (eₑ - p * β * cosθ)
    jac = (p / (2.0 *
            sqrt((1.0 + (β * cosθ)^2) * eₑ^2 - (1.0 + β^2 * (-1.0 + cosθ^2)) * me^2 -
                 2.0 * β * cosθ * eₑ * p) *
            γ))

    dnde = 0.0
    if (fs == "total" || fs == "π⁺ π⁻")
        dnde += 2.0 * pwππ * positron_spectrum_charged_pion(eₑ_srf, mod.ms / 2.0)
    elseif (fs == "total" || fs == "μ⁺ μ⁻")
        dnde += 2.0 * pwμμ * positron_spectrum_muon(eₑ_srf, mod.ms / 2.0)
    end
    jac * dnde
end

"""
    dndeₑ_ss(eₑ::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into scalar mediators.
"""
function dndeₑ_ss(eₑ::Real, es::Real, mod::AbstractScalarMediator; fs::String="total")
    es < ms && return 0.0

    lines_contrib = 0.0
    β = sqrt(1 - (ms / es)^2)

    r = sqrt(1.0 - 4.0 * me^2 / ms^2)
    e₊ = es * (1.0 + r * β) / 2.0
    e₋ = es * (1.0 - r * β) / 2.0
    result = 0.0

    if e₋ <= eₑ <= e₊
        lines_contrib = Γ_s_to_ff(mod, me)  / (es * β)
    end

    fs == "e⁺ e⁻" && return lines_contrib

    if fs == "total" || fs == "π⁺ π⁻" || fs == "μ⁺ μ⁻"
        nodes, weights = gausslegendre(50)
        f(cosθ) = dndeₑ_ss_integrand(cosθ, eₑ, es, mod, fs)
        result = sum(weights .* f.(nodes))
        return result + lines_contrib
    end

    return result
end
