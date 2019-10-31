abstract type AbstractScalarMediator end

mutable struct ScalarMediator <: AbstractScalarMediator
    mχ::Float64
    ms::Float64
    gsχχ::Float64
    gsff::Float64
    gsGG::Float64
    gsFF::Float64
    Λ::Float64

    widths::Float64
    vs::Float64
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
    ScalarMediator(mχ, ms, gsχχ, gsff, gsGG, gsFF, Λ)
end

"""
    list_annihilation_final_states(model::AbstractScalarMediator)

Return the availible final states of the `ScalarMediator` model.
"""
function list_annihilation_final_states(model::AbstractScalarMediator)
    ["μ⁺ μ⁻", "e⁺ e⁻", "γ γ", "π⁰ π⁰", "π⁺ π⁻", "s s"]
end

mutable struct HiggsPortal <: AbstractScalarMediator
    mχ::Float64
    ms::Float64
    gsχχ::Float64
    sinθ::Float64
end

function HiggsPortal(mχ::Float64, ms::Float64, gsχχ::Float64, sinθ::Float64)
    HiggsPortal(mχ, ms, gsχχ, sinθ)
end

getproperty(model::HiggsPortal, :gsff) = model.sinθ
getproperty(model::HiggsPortal, :gsGG) = 3 * model.sinθ
getproperty(model::HiggsPortal, :gsFF) = -5 / 6 * model.sinθ
getproperty(model::HiggsPortal, :Λ) = VH

mutable struct HeavyQuark <: AbstractScalarMediator
    mχ::Float64
    ms::Float64
    gsχχ::Float64
    gsQ::Float64
    mQ::Float64
    QQ::Float64
end

getproperty(model::HeavyQuark, :gsff) = 0.0
getproperty(model::HeavyQuark, :gsGG) = model.gsQ
getproperty(model::HeavyQuark, :gsFF) = 0.0
getproperty(model::HeavyQuark, :Λ) = model.mQ

"""
    list_annihilation_final_states(model::HeavyQuark)

Return the availible final states of the `HeavyQuark` model.
"""
function list_annihilation_final_states(model::HeavyQuark)
    ["γ γ", "π⁰ π⁰", "π⁺ π⁻", "s s"]
end
