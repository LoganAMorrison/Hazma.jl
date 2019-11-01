include("theory.jl")
include("constants.jl")
import Base.getproperty

abstract type AbstractScalarMediator <: AbstractTheoryMediator end

mutable struct ScalarMediator <: AbstractScalarMediator
    mχ::Float64
    ms::Float64
    gsχχ::Float64
    gsff::Float64
    gsGG::Float64
    gsFF::Float64
    Λ::Float64
end

function ScalarMediator(
    mχ::Float64,
    ms::Float64,
    gsχχ::Float64,
    gsff::Float64,
    gsGG::Float64,
    gsFF::Float64,
    Λ::Float64,
)
    @assert mχ >= 0.0
    @assert ms >= 0.0
    ScalarMediator(mχ, ms, gsχχ, gsff, gsGG, gsFF, Λ)
end


# Approximately true for all reasonable parameter values
scalar_vev(::AbstractScalarMediator) = 0.0

function list_decay_final_states(model::AbstractScalarMediator)
    ["e⁺ e⁻", "μ⁺ μ⁻", "γ γ", "π⁰ π⁰", "π⁺ π⁻", "χ χ"]
end

function list_annihilation_final_states(model::AbstractScalarMediator)
    ["e⁺ e⁻", "μ⁺ μ⁻", "γ γ", "π⁰ π⁰", "π⁺ π⁻", "s s"]
end

mutable struct HiggsPortal <: AbstractScalarMediator
    mχ::Float64
    ms::Float64
    gsχχ::Float64
    sinθ::Float64
end

function Base.getproperty(model::HiggsPortal, param::Symbol)
    if param == :gsff
        return model.sinθ
    elseif param == :gsGG
        return 3 * model.sinθ
    elseif param == :gsFF
        return -5 / 6 * model.sinθ
    elseif param == :Λ
        return VH
    else
        return getfield(model, param)
    end
end

mutable struct HeavyQuark <: AbstractScalarMediator
    mχ::Float64
    ms::Float64
    gsχχ::Float64
    gsQQ::Float64
    mQ::Float64
    QQ::Float64
end

function Base.getproperty(model::HeavyQuark, param::Symbol)
    if param == :gsff
        return 0.0
    elseif param == :gsGG
        return model.gsQQ
    elseif param == :gsFF
        return 0.0
    elseif param == :Λ
        return model.mQ
    else
        return getfield(model, param)
    end
end

function list_annihilation_final_states(model::HeavyQuark)
    ["γ γ", "π⁰ π⁰", "π⁺ π⁻", "s s"]
end

## Mediator widths

function Γ_med(mod::AbstractScalarMediator, fs::String)
    (Γ_s_to_γγ(mod) + Γ_s_to_π0π0(mod) + Γ_s_to_ππ(mod) + Γ_s_to_χχ(mod) +
     Γ_s_to_ee(mod) + Γ_s_to_μμ(mod))
end

function Γ_s_to_γγ(mod::AbstractScalarMediator)
    (αem^2 * mod.gsFF^2 * mod.ms^3) / (128.0 * mod.Λ^2 * π^3)
end

function Γ_s_to_π0π0(mod::AbstractScalarMediator)
    (mod.ms < 2mπ⁰) && return 0.0
    vs = scalar_vev(mod)
    return (sqrt(-4 * 2 * mπ⁰^2 + mod.ms^2) *
            ((mod.gsGG * (4 * mπ⁰^2 - 2 * mod.ms^2)) /
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

function Γ_s_to_ff(mod::AbstractScalarMediator, mf::Float64)
    mod.ms < 2mf && return 0.0
    return (-(mod.gsff^2 * mf^2 * (2 * mf - mod.ms) * (2 * mf + mod.ms) *
              sqrt(-4 * mf^2 + mod.ms^2)) / (32.0 * mod.ms^2 * π * VH^2))
end

Γ_s_to_ee(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, me)
Γ_s_to_μμ(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, mμ)

### DM self-annihilation cross sections

function σ_χχ_to_ff(e_cm::Float64, mod::AbstractScalarMediator, mf::Float64)
    (e_cm < 2mf || e_cm < 2 * mod.mχ) && return 0.0
    Γs = Γ_med(mod)
    return ((mod.gsFF^2 * mod.gsχχ^2 * mf^2 * (-2 * mod.mχ + e_cm) *
             (2 * mod.mχ + e_cm) * (-4 * mf^2 + e_cm^2)^1.5) /
            (16π * e_cm^2 * √(e_cm^2 - 4 * mod.mχ^2) * VH^2 *
             (mod.ms^4 - 2 * mod.ms^2 * e_cm^2 + e_cm^4 + mod.ms^2 * Γs^2)))
end

σ_χχ_to_ee(e_cm::Float64, mod::AbstractScalarMediator) =
    σ_χχ_to_ff(e_cm, mod, me)
σ_χχ_to_μμ(e_cm::Float64, mod::AbstractScalarMediator) =
    σ_χχ_to_ff(e_cm, mod, mμ)

function σ_χχ_to_γγ(e_cm::Float64, mod::AbstractScalarMediator)
    (e_cm < 2 * mod.mχ) && return 0.0

    Γs = Γ_med(mod)
    rxs = mod.mχ / e_cm
    return ((αem^2 * mod.gsFF^2 * mod.gsχχ^2 * e_cm^4 * √(1 - 4 * rxs^2)) /
            (64 * mod.Λ^2 * π^3 *
             (mod.ms^4 + e_cm^4 + mod.ms^2 * (-2 * e_cm^2 + Γs^2))))
end

function σ_χχ_to_π⁰π⁰(e_cm::Float64, mod::AbstractScalarMediator)
    (e_cm < 2mπ⁰ || e_cm < 2 * mod.mχ) && return 0.0

    Γs = Γ_med(mod)
    vs = scalar_vev(mod)
    rπ⁰s = mπ⁰ / e_cm
    rxs = mod.mχ / e_cm

    return ((mod.gsχχ^2 * √((-1 + 4 * rπ⁰s^2) * (-1 + 4 * rxs^2)) *
             (162 * mod.gsGG * mod.Λ^3 * e_cm^2 * (-1 + 2 * rπ⁰s^2) * VH^2 +
              B0 * (md + mu) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
              (27 * mod.gsFF^2 * mod.Λ^2 * vs *
               (3 * mod.Λ + 4 * mod.gsGG * vs) - 2 * mod.gsGG * VH^2 *
               (27 * mod.Λ^2 - 30 * mod.gsGG * mod.Λ * vs +
                8 * mod.gsGG^2 * vs^2) +
               mod.gsFF *
               (-81 * mod.Λ^3 * VH + 48 * mod.gsGG^2 * mod.Λ * VH * vs^2)))^2) /
            (209952 * mod.Λ^6 * π * VH^4 * (9 * mod.Λ + 4 * mod.gsGG * vs)^2 *
             (mod.ms^4 + e_cm^4 + mod.ms^2 * (-2 * e_cm^2 + Γs^2))))
end

function σ_χχ_to_ππ(e_cm::Float64, mod::AbstractScalarMediator)
    (e_cm < 2.0 * mπ || e_cm < 2.0 * mod.mχ) && return 0.0

    Γs = Γ_med(mod)
    vs = scalar_vev(mod)
    rpis = mπ / e_cm
    rxs = mod.mχ / e_cm

    return ((mod.gsχχ^2 * √((-1 + 4 * rpis^2) * (-1 + 4 * rxs^2)) *
             (162 * mod.gsGG * mod.Λ^3 * e_cm^2 * (-1 + 2 * rpis^2) * VH^2 +
              B0 * (md + mu) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
              (27 * mod.gsFF^2 * mod.Λ^2 * vs *
               (3 * mod.Λ + 4 * mod.gsGG * vs) - 2 * mod.gsGG * VH^2 *
               (27 * mod.Λ^2 - 30 * mod.gsGG * mod.Λ * vs +
                8 * mod.gsGG^2 * vs^2) +
               mod.gsFF *
               (-81 * mod.Λ^3 * VH + 48 * mod.gsGG^2 * mod.Λ * VH * vs^2)))^2) /
            (104976.0 * mod.Λ^6 * π * VH^4 * (9 * mod.Λ + 4 * mod.gsGG * vs)^2 *
             (mod.ms^4 + e_cm^4 + mod.ms^2 * (-2 * e_cm^2 + Γs^2))))
end

function σ_χχ_to_ss(e_cm::Float64, mod::AbstractScalarMediator)
    (e_cm < 2.0 * mod.ms || e_cm < 2.0 * mod.mχ) && return 0.0

    return -((mod.gsχχ^4 * √(-4 * mod.ms^2 + e_cm^2) *
              √(-4 * mod.mχ^2 + e_cm^2) *
              (-2 / (4 * mod.mχ^2 - e_cm^2) -
               (mod.ms^2 - 4 * mod.mχ^2)^2 / ((4 * mod.mχ^2 - e_cm^2) *
                (mod.ms^4 - 4 * mod.ms^2 * mod.mχ^2 + mod.mχ^2 * e_cm^2)) -
               (2 *
                (6 * mod.ms^4 - 32 * mod.mχ^4 + 16 * mod.mχ^2 * e_cm^2 +
                 e_cm^4 - 4 * mod.ms^2 * (4 * mod.mχ^2 + e_cm^2)) *
                atanh((√(-4 * mod.ms^2 + e_cm^2) * √(-4 * mod.mχ^2 + e_cm^2)) /
                      (-2 * mod.ms^2 + e_cm^2))) /
               (√(-4 * mod.ms^2 + e_cm^2) * (-2 * mod.ms^2 + e_cm^2) *
                (-4 * mod.mχ^2 + e_cm^2)^1.5))) / (16.0 * π * e_cm^2))
end

function σ_χχ(e_cm::Float64, mod::AbstractScalarMediator, fs::String)
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
    end
end


# ----------- #
# FSR Spectra #
# ----------- #

function dnde_χχ_to_ffγ(
    eγ::Real,
    e_cm::Real,
    mod::AbstractScalarMediator,
    mf::Real,
)
    e = eγ / e_cm
    rf = mf / e_cm
    s = e_cm^2 - 2 * e_cm * eγ

    (e_cm < 2mf || e_cm < 2 * mod.mχ || s < 4 * mf^2 || e_cm^2 < s) &&
    return 0.0

    return (αem * (2 * (-1 + 4 * rf^2) *
             sqrt((-1 + 2 * e) * (-1 + 2 * e + 4 * rf^2)) +
             4 * (1 + 2 * (-1 + e) * e - 6 * rf^2 + 8 * e * rf^2 + 8 * rf^4) *
             atanh(sqrt(1 + (4 * rf^2) / (-1 + 2 * e))))) /
           (e * (1 - 4 * rf^2)^1.5 * π * e_cm)
end

function dnde_χχ_to_ππγ(eγ::Real, e_cm::Real, mod::AbstractScalarMediator)
    μπ = mπ / e_cm
    x = 2eγ / e_cm

    (x < 0 || 1 - 4 * μπ^2 < x || 2 * mod.mχ > e_cm) && return 0.0

    dynamic = (-2 * sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x) +
               (-1 + 2 * μπ^2 + x) *
               log((1 - x - sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x))^2 /
                   (-1 + x - sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x))^2)) / x

    coeff = αem / (sqrt(1 - 4 * μπ^2) * π)

    2.0 * dynamic * coeff / e_cm
end

# ---------------- #
# Positron Spectra #
# ---------------- #

"""
    dnde_pos_ππ(ep::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into charged pions.
"""
dnde_pos_ππ(ep::Real, e_cm::Real, ::AbstractScalarMediator) =
    positron_spectrum_charged_pion(ep, e_cm / 2.0)

"""
    dnde_pos_μμ(ep::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into muons.
"""
dnde_pos_μμ(ep::Real, e_cm::Real, ::AbstractScalarMediator) =
    positron_spectrum_muon(ep, e_cm / 2.0)


function dnde_pos_ss_integrand(
    cosθ::Real,
    ep::Real,
    es::Real,
    mod::AbstractScalarMediator,
    fs::String,
)
    ep < me && return 0.0

    pwμμ = Γ_s_to_ff(mod, mμ)
    pwππ = Γ_s_to_ππ(mod)

    p = sqrt(ep^2 - me^2)
    γ = es / ms
    β = sqrt(1 - (ms / es)^2)
    epsrf = γ * (ep - p * β * cosθ)
    jac = (p / (2 *
            sqrt((1 + (β * cosθ)^2) * ep^2 - (1 + β^2 * (-1 + cosθ^2)) * me^2 -
                 2 * β * cosθ * ep * p) *
            γ))

    dnde = 0.0
    if (fs == "total" || fs == "π⁺ π⁻")
        dnde += 2 * pwππ * positron_spectrum_charged_pion(epsrf, mod.ms / 2)
    elseif (fs == "total" || fs == "μ⁺ μ⁻")
        dnde += 2 * pwμμ * positron_spectrum_muon(epsrf, mod.ms / 2)
    end
    jac * dnde
end

"""
    dnde_pos_ss(ep::Real, e_cm::Real)

Positron/electron spectrum from dark matter annihilating into scalar mediators.
"""
function dnde_pos_ss(ep::Real, es::Real, mod::AbstractScalarMediator; fs::String="total")
    es < ms && return 0.0

    lines_contrib = 0.0
    β = sqrt(1 - (ms / es)^2)

    r = sqrt(1 - 4 * me^2 / ms^2)
    eplus = es * (1 + r * β) / 2
    eminus = es * (1 - r * β) / 2
    result = 0.0

    if eminus <= ep <= eplus
        lines_contrib = Γ_s_to_ff(mod, me)  / (es * β)
    end

    fs == "e⁺ e⁻" && return lines_contrib

    if fs == "total" || fs == "π⁺ π⁻" || fs == "μ⁺ μ⁻"
        nodes, weights = gausslegendre(50)
        f(cosθ) = dnde_pos_ss_integrand(cosθ, ep, es, mod, fs)
        result = sum(weights .* f.(nodes))
        return result + lines_contrib
    end

    return result
end
