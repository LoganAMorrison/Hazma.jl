# Scalar mediator widths in the scalar mediator mod.

"""
    Γ_med(mod::AbstractScalarMediator, fs::String)

Compute the partial or total decay width of a scalar mediator into the final
state `fs` from model `mod`
"""
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

"""
    Γ_s_to_γγ(mod::AbstractScalarMediator)

Compute the partial decay width of a scalar mediator into photons from model
`mod`.
"""
function Γ_s_to_γγ(mod::AbstractScalarMediator)
    (αem^2 * mod.gsFF^2 * mod.ms^3) / (64 * mod.Λ^2 * π^3)
end

"""
    Γ_s_to_π⁰π⁰(mod::AbstractScalarMediator)

Compute the partial decay width of a scalar mediator into neutral pions from
model `mod`.
"""
function Γ_s_to_π⁰π⁰(mod::AbstractScalarMediator)
    (mod.ms < 2mπ⁰) && return 0.0
    vs = scalar_vev(mod)

    ((sqrt(-4 * mπ⁰^2 + mod.ms^2) *
      (-162 * mod.gsGG * mod.Λ^3 * (-2 * mπ⁰^2 + mod.ms^2) * VH^2 +
       B0 * (md + mu) * (9 * mod.Λ + 4 * mod.gsGG * vs) *
       (-3 * mod.Λ * VH + 3 * mod.gsff * mod.Λ * vs + 2 * mod.gsGG * VH * vs) *
       (2 * mod.gsGG * VH * (9 * mod.Λ - 4 * mod.gsGG * vs) +
        9 * mod.gsff * mod.Λ * (3 * mod.Λ + 4 * mod.gsGG * vs)))^2) /
     (209952 * mod.Λ^6 * mod.ms^2 * π * VH^4 *
      (9 * mod.Λ + 4 * mod.gsGG * vs)^2))
end

"""
    Γ_s_to_ππ(mod::AbstractScalarMediator)

Compute the partial decay width of a scalar mediator into charged pions from
model `mod`.
"""
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

"""
    Γ_s_to_χχ(mod::AbstractScalarMediator)

Compute the partial decay width of a scalar mediator into dark matter from
model `mod`.
"""
function Γ_s_to_χχ(mod::AbstractScalarMediator)
    (mod.ms < 2 * mod.mχ) && return 0.0

    (mod.gsχχ^2 * (mod.ms^2 - 4 * mod.mχ^2)^1.5) / (8 * mod.ms^2 * π)
end

"""
    Γ_s_to_ff(mod::AbstractScalarMediator)

Compute the partial decay width of a scalar mediator into fermions of mass `mf`
from model `mod`. Only the electron or muon should be used in this function.
However, no checks are performed.
"""
function Γ_s_to_ff(mod::AbstractScalarMediator, mf::Real)
    mod.ms < 2mf && return 0.0
    gsll = mod.gsff * (mf / VH)

    (gsll^2 * (-4 * mf^2 + mod.ms^2)^1.5) / (8 * mod.ms^2 * π)
end

"""
    Γ_s_to_ee(mod::AbstractScalarMediator)

Compute the partial decay width of a scalar mediator into electrons from model
`mod`.
"""
Γ_s_to_ee(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, me)

"""
    Γ_s_to_μμ(mod::AbstractScalarMediator)

Compute the partial decay width of a scalar mediator into muons from model
`mod`.
"""
Γ_s_to_μμ(mod::AbstractScalarMediator) = Γ_s_to_ff(mod, mμ)

"""
    Γ_s_to_ee(mod::HeavyQuark)

Compute the partial decay width of a scalar mediator into electrons from HeavyQuark mod. Note this is zero since there are no coupling between the
scalar mediator and leptons in this mod.
"""
Γ_s_to_ee(mod::HeavyQuark) = 0.0;

"""
    Γ_s_to_μμ(mod::HeavyQuark)

Compute the partial decay width of a scalar mediator into muons from HeavyQuark mod. Note this is zero since there are no coupling between the scalar
mediator and leptons in this mod.
"""
Γ_s_to_μμ(mod::HeavyQuark) = 0.0;
