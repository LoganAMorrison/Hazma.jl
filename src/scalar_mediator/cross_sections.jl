# Cross sections in the scalar mediator model.

"""
    σ_χχ(e_cm::Real, mod::AbstractScalarMediator, fs::String)

Compute the cross section of χχ̄ → `fs` in the scalar mediator model `mod`. Use
`fs = "total"` to sum over all final states.
"""
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
    elseif fs == "total"
        return sum(σ_χχ(e_cm, mod, fs) for fs in list_annihilation_final_states(mod))
    else
        return zero(typeof(e_cm))
    end
end

"""
    σ_χχ_to_ll(e_cm::Real, mod::AbstractScalarMediator, ml::Real)

Compute the cross section of χχ̄ → ll̄ in the scalar mediator model `mod` given
a final state lepton mass `ml`.
"""
function σ_χχ_to_ll(e_cm::Real, mod::AbstractScalarMediator, ml::Real)
    (e_cm < 2 * mf || e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))

    (mod.gsff^2 * mod.gsχχ^2 * ml^2 * (-4 * ml^2 + e_cm^2)^1.5 *
     sqrt(-4 * mod.mχ^2 + e_cm^2)) /
    (16 * π * e_cm^2 * VH^2 * ((mod.ms^2 - e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))
end

"""
    σ_χχ_to_ee(e_cm::Real, mod::AbstractScalarMediator)

Compute the cross section of χχ̄ → e⁺e⁻ in the scalar mediator model `mod`.
"""
σ_χχ_to_ee(e_cm::Real, mod::AbstractScalarMediator) = σ_χχ_to_ll(e_cm, mod, me)

"""
    σ_χχ_to_μμ(e_cm::Real, mod::AbstractScalarMediator)

Compute the cross section of χχ̄ → μ⁺μ⁻ in the scalar mediator model `mod`.
"""
σ_χχ_to_μμ(e_cm::Real, mod::AbstractScalarMediator) = σ_χχ_to_ll(e_cm, mod, mμ)

"""
    σ_χχ_to_ee(e_cm::Real, mod::HeavyQuark)

Compute the cross section of χχ̄ → e⁺e⁻ in the HeavyQuark model `mod`. Note
this is zero since there are no couplings between the scalar mediator and
leptons in the HeavyQuark model.
"""
σ_χχ_to_ee(e_cm::Real, mod::HeavyQuark) = zero(e_cm)

"""
    σ_χχ_to_μμ(e_cm::Real, mod::HeavyQuark)

Compute the cross section of χχ̄ → μ⁺μ⁻ in the HeavyQuark model `mod`. Note
this is zero since there are no couplings between the scalar mediator and
leptons in the HeavyQuark model.
"""
σ_χχ_to_μμ(e_cm::Real, mod::HeavyQuark) = zero(e_cm)

"""
    σ_χχ_to_γγ(e_cm::Real, mod::AbstractScalarMediator)

Compute the cross section of χχ̄ → γγ in the scalar mediator model `mod`.
"""
function σ_χχ_to_γγ(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))
    Γs = mod.Γ_med

    (αem^2 * mod.gsFF^2 * mod.gsχχ^2 * e_cm^3 * sqrt(-4 * mod.mχ^2 + e_cm^2)) /
    (128 * mod.Λ^2 * π^3 * ((mod.ms^2 - e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))
end

"""
    σ_χχ_to_π⁰π⁰(e_cm::Real, mod::AbstractScalarMediator)

Compute the cross section of χχ̄ → π⁰π⁰ in the scalar mediator model `mod`.
"""
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

"""
    σ_χχ_to_ππ(e_cm::Real, mod::AbstractScalarMediator)

Compute the cross section of χχ̄ → π⁺π⁻ in the scalar mediator model `mod`.
"""
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

"""
    σ_χχ_to_ss(e_cm::Real, mod::AbstractScalarMediator)

Compute the cross section of χχ̄ → ss in the scalar mediator model `mod`.
"""
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

"""
    σ_χχ_to_χχ(e_cm::Real, mod::AbstractScalarMediator)

Compute the self interaction cross section of χχ̄ → χχ̄ in the scalar mediator
model `mod`.
"""
function σ_χχ_to_χχ(e_cm::Real, mod::AbstractScalarMediator)
    (e_cm < 2 * mod.mχ) && return zero(typeof(e_cm))

    (mod.gsχχ^4 * (2 *
      ((mod.ms^2 - 4 * mod.mχ^2)^2 * (mod.ms^2 - e_cm^2)^2 +
       2 * mod.ms^2 * e_cm^2 * (2 * mod.ms^2 - e_cm^2) * mod.Γ_med^2 -
       mod.ms^4 * mod.Γ_med^4) *
      (-acot(mod.Γ_med / mod.ms) +
       acot((mod.ms * mod.Γ_med) / (mod.ms^2 - 4 * mod.mχ^2 + e_cm^2))) +
      mod.ms * mod.Γ_med *
      (-2 * (4 * mod.mχ^2 - e_cm^2) *
       (mod.ms^4 + 16 * mod.mχ^4 - 4 * mod.mχ^2 * e_cm^2 + 3 * e_cm^4 +
        mod.ms^2 * (-4 * mod.mχ^2 - 3 * e_cm^2 + mod.Γ_med^2)) +
       ((mod.ms - e_cm) * (mod.ms + e_cm) *
        (2 * mod.ms^4 - 3 * mod.ms^2 * (4 * mod.mχ^2 + e_cm^2) +
         4 * mod.mχ^2 * (4 * mod.mχ^2 + e_cm^2)) +
        mod.ms^2 * (2 * mod.ms^2 - 4 * mod.mχ^2 + e_cm^2) * mod.Γ_med^2) *
       log((mod.ms^2 * (mod.ms^2 + mod.Γ_med^2)) /
           ((mod.ms^2 - 4 * mod.mχ^2 + e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))))) /
    (32 * mod.ms * pi * e_cm^2 * (-4 * mod.mχ^2 + e_cm^2) * mod.Γ_med *
     ((mod.ms^2 - e_cm^2)^2 + mod.ms^2 * mod.Γ_med^2))
end

"""
    thermal_cross_section(x::Real, mod::AbstractScalarMediator)

Compute the thermally average cross section for the dark matter particle of the
given model `mod`.
"""
function thermal_cross_section(x::Real, mod::AbstractScalarMediator)
    x > 300 && return 0.0

    integrand(z) =
        (z^2 * (z^2 - 4) * besselk1(x * z) * σ_χχ(mod.mχ * z, mod, "total"))

    #TODO: add breakpoints at: ms / mχ and 2ms / mχ
    x / (2 * besselk2(x))^2 * quadgk(integrand, 2.0, max(50.0 / x, 150))[1]
end
