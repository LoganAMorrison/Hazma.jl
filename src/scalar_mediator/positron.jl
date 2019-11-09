# Positron spectra in the scalar mediator model.

"""
    dndeₑ_χχ(eₑ::Real, e_cm::Real, mod::AbstractScalarMediator, fs::String)

Compute the electron/positron spectrum from χχ̄ → `fs` given a center of mass
energy `e_cm` and scalar mediator model `mod`.
"""
function dndeₑ_χχ(eₑ::Real, e_cm::Real, mod::AbstractScalarMediator, fs::String)
    if fs == "π⁺ π⁻"
        return dndeₑ_χχ_to_ππ(eₑ, e_cm, mod)
    elseif fs == "μ μ"
        return dndeₑ_χχ_to_μμ(eₑ, e_cm, mod)
    elseif fs == "s s"
        return dndeₑ_χχ_to_ss(eₑ, e_cm, mod)
    else
        return zero(typeof(eₑ))
    end
end

"""
    lines_e(e_cm::Real, mod::AbstractScalarMediator)

Compute the electron/positron lines from χχ̄ → e⁺e⁻ given the center of mass
energy `e_cm` and scalar mediator model `mod`. Returns a dictionary of the
location of the lines.
"""
function lines_e(e_cm::Real, mod::AbstractScalarMediator)
    Dict{String,T}("e⁺ e⁻" => e_cm / 2)
end

"""
    dndeₑ_χχ_to_ππ(eₑ::Real, e_cm::Real)

Compute the electron/positron spectrum from χχ̄ → π⁺π⁻→ X + e⁺/e⁻ from the
decay of a charged pion given a electron/positron energy `eₑ` and center of mass
energy `e_cm`.
"""
dndeₑ_χχ_to_ππ(eₑ::Real, e_cm::Real, ::AbstractScalarMediator) =
    dndeₑ_π_decay(eₑ, e_cm / 2)

"""
    dndeₑ_χχ_to_μμ(eₑ::Real, e_cm::Real)

Compute the electron/positron spectrum from χχ̄ → μ⁺μ⁻→ X + e⁺/e⁻ from the
decay of a muon given a electron/positron energy `eₑ` and center of mass
energy `e_cm`.
"""
dndeₑ_χχ_to_μμ(eₑ::Real, e_cm::Real, ::AbstractScalarMediator) =
    dndeₑ_μ_decay(eₑ, e_cm / 2)


"""
    dndeₑ_χχ_to_ss(eₑ::Real, e_cm::Real)

Compute the electron/positron spectrum from χχ̄ → ss → X + e⁺/e⁻ from the
decay of a scalar mediator (into μ⁺μ⁻, π⁺π⁻ or e⁺e⁻) given a electron/positron
energy `eₑ` and center of mass energy `e_cm`.
"""
function dndeₑ_χχ_to_ss(
    eₑ::Real,
    e_cm::Real,
    mod::AbstractScalarMediator;
    fs::String = "total",
)
    (e_cm < 2 * mod.mχ || e_cm < 2 * mod.ms) && return zero(typeof(eₑ))
    mod.Γ_med == 0 && return zero(typeof(eₑ))

    es = e_cm / 2

    # line contribution from s → e⁺ e⁻: factor of 2 for 2 electrons per S
    Γee = Γ_s_to_ee(mod) / mod.Γ_med
    lines_contrib = 2 * Γee *
                    boosted_delta_function(es, mod.ms, eₑ, me, mod.ms / 2)

    fs == "e⁺ e⁻" && return lines_contrib
    Γππ = Γ_s_to_ππ(mod) / mod.Γ_med
    Γμμ = Γ_s_to_μμ(mod) / mod.Γ_med
    # factors of 2 for two π (μ) per S
    dndeₑ_π_decay_srf(eₑ_srf) = 2 * Γππ * dndeₑ_π_decay(eₑ_srf, mod.ms / 2)
    dndeₑ_μ_decay_srf(eₑ_srf) = 2 * Γμμ * dndeₑ_μ_decay(eₑ_srf, mod.ms / 2)
    spectrum_rf(eₑ_srf) = dndeₑ_π_decay_srf(eₑ_srf) + dndeₑ_μ_decay_srf(eₑ_srf)

    if fs == "total"
        result = boost_spectrum(spectrum_rf, es, mod.ms, eₑ, me)
    elseif fs == "π⁺ π⁻"
        result = boost_spectrum(dndeₑ_π_decay_srf, es, mod.ms, eₑ, me)
    elseif fs == "μ⁺ μ⁻"
        result = boost_spectrum(dndeₑ_π_decay_srf, es, mod.ms, eₑ, me)
    else
        result = zero(typeof(eₑ))
    end

    # factor of 2 for 2 final state S's
    2 * (result + lines_contrib)
end
