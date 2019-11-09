# Model definitions for various scalar mediator models.

import Base.getproperty

abstract type AbstractScalarMediator <: AbstractTheoryMediator end

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
