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
    # TODO: add inner constructor to validate args
end

# Approximately true for all reasonable parameter values
scalar_vev(AbstractScalarMediator) = 0.0

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

# TODO: implement partial widths!!!
function Γ_med(mod::AbstractScalarMediator, fs::String)
    return 0.0
end

### DM self-annihilation cross sections

function σ_χχ_to_ff(e_cm::Float64, mod::AbstractScalarMediator, mf::Float64)
    if e_cm < 2mf || e_cm < 2mod.mχ
        return 0.
    else
        Γs = Γ_med(mod)
        return ((mod.gsff^2 * mod.gsχχ^2 * mf^2 * (-2mod.mχ + e_cm) * (2mod.mχ + e_cm)
                * (-4mf^2 + e_cm^2)^1.5) / (16π * e_cm^2 * √(e_cm^2 - 4mod.mχ^2) * VH^2
                * (mod.ms^4 - 2mod.ms^2 * e_cm^2 + e_cm^4 + mod.ms^2 * Γs^2)))
    end
end

σ_χχ_to_ee(e_cm::Float64, mod::AbstractScalarMediator) = σ_χχ_to_ff(e_cm, mod, me)
σ_χχ_to_μμ(e_cm::Float64, mod::AbstractScalarMediator) = σ_χχ_to_ff(e_cm, mod, mμ)

function σ_χχ_to_γγ(e_cm::Float64, mod::AbstractScalarMediator)
    if e_cm < 2mod.mχ
        return 0.0
    else
        Γs = Γ_med(mod)
        rxs = mod.mχ / e_cm
        return ((αem^2 * mod.gsFF^2 * mod.gsχχ^2 * e_cm^4 * √(1 - 4rxs^2))
                / (64mod.Λ^2 * π^3 * (mod.ms^4 + e_cm^4 + mod.ms^2 * (-2e_cm^2 + Γs^2))))
    end
end

function σ_χχ_to_π⁰π⁰(e_cm::Float64, mod::AbstractScalarMediator)
    if e_cm < 2mπ⁰ || e_cm < 2mod.mχ
        return 0.
    else
        Γs = Γ_med(mod)
        vs = scalar_vev(mod)
        rπ⁰s = mπ⁰ / e_cm
        rxs = mod.mχ / e_cm

        return ((mod.gsχχ^2 * √((-1 + 4rπ⁰s^2) * (-1 + 4rxs^2)) *
                 (162mod.gsGG * mod.Λ^3 * e_cm^2 *
                  (-1 + 2rπ⁰s^2) * VH^2 + b0 * (mdq + muq) *
                  (9mod.Λ + 4mod.gsGG * vs) *
                  (27mod.gsff^2 * mod.Λ^2 * vs * (3mod.Λ + 4mod.gsGG * vs) -
                   2 * mod.gsGG * VH^2 *
                   (27mod.Λ^2 - 30 * mod.gsGG * mod.Λ * vs +
                    8 * mod.gsGG^2 * vs^2) + mod.gsff *
                   (-81 * mod.Λ^3 * VH +
                    48 * mod.gsGG^2 * mod.Λ * VH * vs^2)))^2) /
                (209952 * mod.Λ^6 * π * VH^4 *
                 (9 * mod.Λ + 4 * mod.gsGG * vs)^2 *
                 (mod.ms^4 + e_cm^4 + mod.ms^2 *
                  (-2 * e_cm^2 + Γs^2))))
    end
end

function σ_χχ_to_ππ(e_cm::Float64, mod::AbstractScalarMediator)
    if e_cm < 2.0 * mπ || e_cm < 2.0 * mod.mχ
        return 0.0
    else
        Γs = Γ_med(mod)
        vs = scalar_vev(mod)
        rpis = mπ / e_cm
        rxs = mod.mχ / e_cm

        return ((mod.gsχχ^2 * √((-1 + 4 * rpis^2) *
                                   (-1 + 4 * rxs^2)) *
                 (162 * mod.gsGG * mod.Λ^3 * e_cm^2 *
                  (-1 + 2 * rpis^2) * VH^2 +
                  b0 * (mdq + muq) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
                  (27 * mod.gsff^2 * mod.Λ^2 * vs *
                   (3 * mod.Λ + 4 * mod.gsGG * vs) -
                   2 * mod.gsGG * VH^2 *
                   (27 * mod.Λ^2 - 30 * mod.gsGG * mod.Λ * vs +
                    8 * mod.gsGG^2 * vs^2) + mod.gsff *
                   (-81 * mod.Λ^3 * VH +
                    48 * mod.gsGG^2 * mod.Λ * VH * vs^2)))^2) /
                (104976.0 * mod.Λ^6 * π * VH^4 *
                 (9 * mod.Λ + 4 * mod.gsGG * vs)^2 *
                 (mod.ms^4 + e_cm^4 + mod.ms^2 *
                  (-2 * e_cm^2 + Γs^2))))
    end
end

function σ_χχ_to_ss(e_cm::Float64, mod::AbstractScalarMediator)
    if e_cm < 2.0 * mod.ms || e_cm < 2.0 * mod.mχ
        return 0.0
    else
        return -((mod.gsχχ^4 * √(-4 * mod.ms^2 + e_cm^2) *
                √(-4 * mod.mχ^2 + e_cm^2) *
                (-2 / (4 * mod.mχ^2 - e_cm^2) -
                (mod.ms^2 - 4 * mod.mχ^2)^2 /
                ((4 * mod.mχ^2 - e_cm^2) *
                (mod.ms^4 - 4 * mod.ms^2 * mod.mχ^2 + mod.mχ^2 * e_cm^2)) -
                (2 * (6 * mod.ms^4 - 32 * mod.mχ^4 + 16 * mod.mχ^2 * e_cm^2 +
                e_cm^4 - 4 * mod.ms^2 *
                (4 * mod.mχ^2 + e_cm^2)) *
                atanh((√(-4 * mod.ms^2 + e_cm^2) *
                √(-4 * mod.mχ^2 + e_cm^2)) /
                    (-2 * mod.ms^2 + e_cm^2))) /
                (√(-4 * mod.ms^2 + e_cm^2) *
                (-2 * mod.ms^2 + e_cm^2) *
                (-4 * mod.mχ^2 + e_cm^2)^1.5))) /
                (16.0 * π * e_cm^2))
    end
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
