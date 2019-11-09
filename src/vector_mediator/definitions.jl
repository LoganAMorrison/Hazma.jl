# Model definitions for various scalar mediator models.

import Base.getproperty

abstract type AbstractVectorMediator <: AbstractTheoryMediator end

@make_mediator_theory VectorMediator [
    (:mχ, nonnegative),
    (:mv, nonnegative),
    :gvχχ,
    :gvuu,
    :gvdd,
    :gvss,
    :gvee,
    :gvμμ,
] AbstractVectorMediator

@make_mediator_theory KineticMixing [
    (:mχ, nonnegative),
    (:mv, nonnegative),
    :gvχχ,
    :ϵ,
] AbstractVectorMediator


function Base.getproperty(mod::KineticMixing, param::Symbol)
    if param == :gvuu
        return -mod.ϵ * sqrt(4 * π * αem) * (2 / 3)
    elseif param == :gvdd
        return -mod.ϵ * sqrt(4 * π * αem) * (-1 / 3)
    elseif param == :gvss
        return -mod.ϵ * sqrt(4 * π * αem) * (-1 / 3)
    elseif param == :gvee
        return -mod.ϵ * sqrt(4 * π * αem) * (-1)
    elseif param == :gvμμ
        return -mod.ϵ * sqrt(4 * π * αem) * (-1)
    else
        return getfield(mod, param)
    end
end

function list_decay_final_states(::AbstractVectorMediator)
    ["e⁺ e⁻", "μ⁺ μ⁻", "π⁰ γ", "π⁺ π⁻", "χ χ"]
end

function list_annihilation_final_states(::AbstractVectorMediator)
    ["e⁺ e⁻", "μ⁺ μ⁻", "π⁰ γ", "π⁰ v", "π⁺ π⁻", "v v"]
end
