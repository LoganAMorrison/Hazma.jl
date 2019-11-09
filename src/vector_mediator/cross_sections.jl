"""
    σ_χχ(e_cm::Real, mod::AbstractVectorMediator, fs::String)

Compute the cross section of χχ̄ → `fs` in the vector mediator model `mod`. Use
`fs = "all"` to return a dictionary of all of the final state cross sections,
barring the self interaction cross section.
"""
function σ_χχ(e_cm::Real, mod::AbstractVectorMediator, fs::String)
    if fs == "e⁺ e⁻"
        return σ_χχ_to_ee(e_cm, mod)
    elseif fs == "μ⁺ μ⁻"
        return σ_χχ_to_μμ(e_cm, mod)
    elseif fs == "π⁰ γ"
        return σ_χχ_to_π⁰γ(e_cm, mod)
    elseif fs == "π⁰ v"
        return σ_χχ_to_π⁰v(e_cm, mod)
    elseif fs == "π⁺ π⁻"
        return σ_χχ_to_ππ(e_cm, mod)
    elseif fs == "v v"
        return σ_χχ_to_vv(e_cm, mod)
    elseif fs == "χ χ"
        return σ_χχ_to_χχ(e_cm, mod)
    elseif fs == "all"
        σs = Dict(fs => σ_χχ(e_cm, mod, fs) for fs in list_annihilation_final_states(mod))
        σs["total"] = sum(values(σs))
        return σs
    elseif fs == "total"
        return sum(σ_χχ(e_cm, mod, fs) for fs in list_annihilation_final_states(mod))
    else
        return zero(typeof(e_cm))
    end
end

"""
    σ_χχ_to_ee(e_cm::Real, mod::AbstractVectorMediator)

Compute the cross section of χχ̄ → e⁺e⁻ in the vector mediator model `mod`.
"""
function σ_χχ_to_ee(e_cm::Real, mod::AbstractVectorMediator)
    (e_cm < 2 * mod.mχ || e_cm < 2 * me) && return zero(e_cm)

    (mod.gvee^2 * mod.gvχχ^2 * sqrt(-4 * me^2 + e_cm^2) * (2 * me^2 + e_cm^2) *
     (2 * mod.mχ^2 + e_cm^2)) /
    (12 * π * e_cm^2 * sqrt(-4 * mod.mχ^2 + e_cm^2) *
     ((mod.mv^2 - e_cm^2)^2 + mod.mv^2 * mod.Γ_med^2))
end

"""
    σ_χχ_to_μμ(e_cm::Real, mod::AbstractVectorMediator)

Compute the cross section of χχ̄ → μ⁺μ⁻ in the vector mediator model `mod`.
"""
function σ_χχ_to_μμ(e_cm::Real, mod::AbstractVectorMediator)
    (e_cm < 2 * mod.mχ || e_cm < 2 * me) && return zero(e_cm)

    (mod.gvμμ^2 * mod.gvχχ^2 * sqrt(-4 * mμ^2 + e_cm^2) * (2 * mμ^2 + e_cm^2) *
     (2 * mod.mχ^2 + e_cm^2)) /
    (12 * π * e_cm^2 * sqrt(-4 * mod.mχ^2 + e_cm^2) *
     ((mod.mv^2 - e_cm^2)^2 + mod.mv^2 * mod.Γ_med^2))
end

"""
    σ_χχ_to_ππ(e_cm::Real, mod::AbstractVectorMediator)

Compute the cross section of χχ̄ → π⁺π⁻ in the vector mediator model `mod`.
"""
function σ_χχ_to_ππ(e_cm::Real, mod::AbstractVectorMediator)
    (e_cm < 2 * mod.mχ || e_cm < 2 * mπ) && return zero(e_cm)

    ((mod.gvdd - mod.gvuu)^2 * mod.gvχχ^2 * (-4 * mπ^2 + e_cm^2)^1.5 *
     (2 * mod.mχ^2 + e_cm^2)) /
    (48 * π * e_cm^2 * sqrt(-4 * mod.mχ^2 + e_cm^2) *
     ((mod.mv^2 - e_cm^2)^2 + mod.mv^2 * mod.Γ_med^2))
end

"""
    σ_χχ_to_π⁰γ(e_cm::Real, mod::AbstractVectorMediator)

Compute the cross section of χχ̄ → π⁰γ in the vector mediator model `mod`.
"""
function σ_χχ_to_π⁰γ(e_cm::Real, mod::AbstractVectorMediator)
    (e_cm < 2 * mod.mχ || e_cm < mπ⁰) && return zero(e_cm)

    (αem * (mod.gvdd + 2 * mod.gvuu)^2 * mod.gvχχ^2 * ((mπ⁰^2 - e_cm^2)^2)^1.5 *
     (2 * mod.mχ^2 + e_cm^2)) /
    (3456 * fπ^2 * π^4 * e_cm^3 * sqrt(-4 * mod.mχ^2 + e_cm^2) *
     ((mod.mv^2 - e_cm^2)^2 + mod.mv^2 * mod.Γ_med^2))
end

"""
    σ_χχ_to_π⁰v(e_cm::Real, mod::AbstractVectorMediator)

Compute the cross section of χχ̄ → π⁰v in the vector mediator model `mod`.
"""
function σ_χχ_to_π⁰v(e_cm::Real, mod::AbstractVectorMediator)
    (e_cm < 2 * mod.mχ || e_cm < mπ⁰ + mod.mv) && return zero(e_cm)

    ((mod.gvdd - mod.gvuu)^2 * (mod.gvdd + mod.gvuu)^2 * mod.gvχχ^2 *
     ((mπ⁰ - mod.mv - e_cm) * (mπ⁰ + mod.mv - e_cm) * (mπ⁰ - mod.mv + e_cm) *
      (mπ⁰ + mod.mv + e_cm))^1.5 *
     (2 * mod.mχ^2 + e_cm^2)) /
    (1536 * fπ^2 * π^5 * e_cm^3 * sqrt(-4 * mod.mχ^2 + e_cm^2) *
     ((mod.mv^2 - e_cm^2)^2 + mod.mv^2 * mod.Γ_med^2))
end

"""
    σ_χχ_to_vv(e_cm::Real, mod::AbstractVectorMediator)

Compute the cross section of χχ̄ → vv in the vector mediator model `mod`.
"""
function σ_χχ_to_vv(e_cm::Real, mod::AbstractVectorMediator)
    (e_cm < 2 * mod.mχ || e_cm < 2 * mod.mv) && return zero(e_cm)

    (mod.gvχχ^4 *
     ((2 * sqrt(-4 * mod.mv^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2) *
       (2 * mod.mv^4 + 4 * mod.mχ^4 + mod.mχ^2 * e_cm^2)) /
      (mod.mv^4 - 4 * mod.mv^2 * mod.mχ^2 + mod.mχ^2 * e_cm^2) +
      ((4 * mod.mv^4 - 8 * mod.mv^2 * mod.mχ^2 - 8 * mod.mχ^4 +
        4 * mod.mχ^2 * e_cm^2 + e_cm^4) * log((-2 * mod.mv^2 + e_cm^2 +
            sqrt(-4 * mod.mv^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2))^2 /
           (2 * mod.mv^2 - e_cm^2 +
            sqrt(-4 * mod.mv^2 + e_cm^2) * sqrt(-4 * mod.mχ^2 + e_cm^2))^2)) /
      (2 * mod.mv^2 - e_cm^2))) / (64 * mod.mχ^2 * π * e_cm^2 - 16 * π * e_cm^4)
end

"""
    σ_χχ_to_χχ(e_cm::Real, mod::AbstractVectorMediator)

Compute the self interaction cross section of χχ̄ → χχ̄ in the vector mediator
model `mod`.
"""
function σ_χχ_to_χχ(e_cm::Real, mod::AbstractVectorMediator)
    (e_cm < 2 * mod.mχ) && return zero(e_cm)

    (mod.gvχχ^4 * (((4 * mod.mχ^2 - e_cm^2) * (3 * mod.mv^4 -
        2 * (4 * mod.mχ^4 + 10 * mod.mχ^2 * e_cm^2 + 7 * e_cm^4) +
        3 * mod.mv^2 *
        (4 * mod.mχ^2 + 3 * (e_cm - mod.Γ_med) * (e_cm + mod.Γ_med)))) / 3.0 +
      (((mod.mv - e_cm)^2 * (mod.mv + e_cm)^2 *
        (mod.mv^4 + 2 * mod.mv^2 * e_cm^2 + 2 * (-2 * mod.mχ^2 + e_cm^2)^2) +
        mod.mv^2 *
        (6 * mod.mv^4 + 8 * (mod.mv - mod.mχ) * (mod.mv + mod.mχ) * e_cm^2 -
         e_cm^4) *
        mod.Γ_med^2 - 3 * mod.mv^4 * mod.Γ_med^4) *
       (-acot(mod.Γ_med / mod.mv) +
        acot((mod.mv * mod.Γ_med) / (mod.mv^2 - 4 * mod.mχ^2 + e_cm^2))) +
       2 * mod.mv * mod.Γ_med *
       (-2 * mod.mχ^4 * e_cm^2 + e_cm^6 -
        mod.mv^4 * (e_cm^2 - 2 * mod.Γ_med^2) +
        mod.mv^2 * (2 * mod.mχ^4 + e_cm^2 * mod.Γ_med^2)) *
       log((mod.mv^2 * (mod.mv^2 + mod.Γ_med^2)) /
           ((mod.mv^2 - 4 * mod.mχ^2 + e_cm^2)^2 + mod.mv^2 * mod.Γ_med^2))) /
      (mod.mv * mod.Γ_med))) / (8 * π * e_cm^2 * (-4 * mod.mχ^2 + e_cm^2) *
     ((mod.mv^2 - e_cm^2)^2 + mod.mv^2 * mod.Γ_med^2))
end
