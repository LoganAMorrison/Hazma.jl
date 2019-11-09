"""
    Γ_med(mod::AbstractVectorMediator, fs::String)

Compute the partial or total decay width of a vector mediator into the final
state `fs` from model `mod`
"""
function Γ_med(mod::AbstractVectorMediator, fs::String)
    if fs == "e⁺ e⁻"
        return Γ_v_to_ee(mod)
    elseif fs == "μ⁺ μ⁻"
        return Γ_v_to_μμ(mod)
    elseif fs == "π⁺ π⁻"
        return Γ_v_to_ππ(mod)
    elseif fs == "χ χ"
        return Γ_v_to_χχ(mod)
    elseif fs == "π⁰ γ"
        return Γ_v_to_π⁰γ(mod)
    elseif fs == "all"
        Γs = Dict(fs => Γ_med(mod, fs) for fs in list_decay_final_states(mod))
        Γs["total"] = sum(values(Γs))
        return Γs
    elseif fs == "total"
        return sum(Γ_med(mod, fs) for fs in list_decay_final_states(mod))
    else
        return 0.0
    end
end

"""
    Γ_s_to_ee(mod::AbstractVectorMediator)

Compute the partial decay width of a vector mediator into electron from model
`mod`.
"""
function Γ_v_to_ee(mod)
    (mod.mv < 2 * me) && return 0.0

    (mod.gvee^2 * sqrt(-4 * me^2 + mod.mv^2) * (2 * me^2 + mod.mv^2)) /
    (12 * mod.mv^2 * π)
end

"""
    Γ_s_to_μμ(mod::AbstractVectorMediator)

Compute the partial decay width of a vector mediator into muons from model
`mod`.
"""
function Γ_v_to_μμ(mod)
    (mod.mv < 2 * mμ) && return 0.0

    (mod.gvμμ^2 * sqrt(-4 * mμ^2 + mod.mv^2) * (2 * mμ^2 + mod.mv^2)) /
    (12 * mod.mv^2 * π)
end

"""
    Γ_s_to_ππ(mod::AbstractVectorMediator)

Compute the partial decay width of a vector mediator into charged pions from
model `mod`.
"""
function Γ_v_to_ππ(mod)
    (mod.mv < 2 * mπ) && return 0.0

    ((mod.gvdd - mod.gvuu)^2 * (-4 * mπ^2 + mod.mv^2)^1.5) / (48 * mod.mv^2 * π)
end

"""
    Γ_s_to_π⁰γ(mod::AbstractVectorMediator)

Compute the partial decay width of a vector mediator into a neutral pion and a
photon from model `mod`.
"""
function Γ_v_to_π⁰γ(mod)
    (mod.mv < mπ) && return 0.0

    (αem * (mod.gvdd + 2 * mod.gvuu)^2 * ((mπ⁰^2 - mod.mv^2)^2)^1.5) /
    (3456 * fπ^2 * mod.mv^3 * π^4)
end

"""
    Γ_s_to_χχ(mod::AbstractVectorMediator)

Compute the partial decay width of a vector mediator into dark matter from
model `mod`.
"""
function Γ_v_to_χχ(mod)
    (mod.mv < 2 * mod.mχ) && return 0.0

    (mod.gvχχ^2 * sqrt(mod.mv^2 - 4 * mod.mχ^2) * (mod.mv^2 + 2 * mod.mχ^2)) /
    (12 * mod.mv^2 * π)
end
