"""
    dndeᵧ_χχ_to_ffγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator, mf::Real)

Compute the photon spectrum from χχ̄ → ff̄γ given a photon energy `eᵧ`, a center
of mass energy `e_cm`, the scalar mediator model `mod` and a final state fermion
mass `mf`.
"""
function dndeᵧ_χχ_to_ffγ(
    eᵧ::Real,
    e_cm::Real,
    mod::AbstractScalarMediator,
    mf::Real,
)
    e = eᵧ / e_cm
    rf = mf / e_cm
    s = e_cm^2 - 2 * e_cm * eᵧ

    (e_cm < 2 * mf || e_cm < 2 * mod.mχ || s < 4 * mf^2 || e_cm^2 < s) &&
    return zero(typeof(eᵧ))
    return (αem * (2 * (-1 + 4 * rf^2) *
             sqrt((-1 + 2 * e) * (-1 + 2 * e + 4 * rf^2)) +
             4 * (1 + 2 * (-1 + e) * e - 6 * rf^2 + 8 * e * rf^2 + 8 * rf^4) *
             atanh(sqrt(1 + (4 * rf^2) / (-1 + 2 * e))))) /
           (e * (1 - 4 * rf^2)^1.5 * π * e_cm)
end

"""
    dndeᵧ_χχ_to_eeγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)

Compute the photon spectrum from χχ̄ → e⁺e⁻γ given a photon energy `eᵧ`, a
center of mass energy `e_cm`, the scalar mediator model `mod`.
"""
dndeᵧ_χχ_to_eeγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator) =
    dndeᵧ_χχ_to_ffγ(eᵧ, e_cm, mod, me)

"""
    dndeᵧ_χχ_to_μμγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)

Compute the photon spectrum from χχ̄ → μ⁺μ⁻γ given a photon energy `eᵧ`, a
center of mass energy `e_cm`, the scalar mediator model `mod`.
"""
dndeᵧ_χχ_to_μμγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator) =
    dndeᵧ_χχ_to_ffγ(eᵧ, e_cm, mod, mμ)

"""
    dndeᵧ_χχ_to_eeγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)

Compute the photon spectrum from χχ̄ → e⁺e⁻γ given a photon energy `eᵧ`, a
center of mass energy `e_cm`, the HeavyQuark model `mod`. Note this is zero
since there are no coupling between the scalar mediator and leptons in the
HeavyQuark model.
"""
dndeᵧ_χχ_to_eeγ(eᵧ::Real, e_cm::Real, mod::HeavyQuark) = zero(eᵧ)

"""
    dndeᵧ_χχ_to_μμγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)

Compute the photon spectrum from χχ̄ → μ⁺μ⁻γ given a photon energy `eᵧ`, a
center of mass energy `e_cm`, the HeavyQuark model `mod`. Note this is zero
since there are no coupling between the scalar mediator and leptons in the
HeavyQuark model.
"""
dndeᵧ_χχ_to_μμγ(eᵧ::Real, e_cm::Real, mod::HeavyQuark) = zero(eᵧ)

"""
    dndeᵧ_χχ_to_ππγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)

Compute the photon spectrum from χχ̄ → π⁺π⁻γ given a photon energy `eᵧ`, a
center of mass energy `e_cm`, the scalar mediator model `mod`.
"""
function dndeᵧ_χχ_to_ππγ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)
    μπ = mπ / e_cm
    x = 2 * eᵧ / e_cm

    (x < 0 || 1 - 4 * μπ^2 < x || 2 * mod.mχ > e_cm) && return zero(typeof(eᵧ))
    dynamic = (-2 * sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x) +
               (-1 + 2 * μπ^2 + x) *
               log((1 - x - sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x))^2 /
                   (-1 + x - sqrt(1 - x) * sqrt(1 - 4 * μπ^2 - x))^2)) / x

    coeff = αem / (sqrt(1 - 4 * μπ^2) * π)

    2 * dynamic * coeff / e_cm
end

"""
    dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)

Compute the photon spectrum from χχ̄ → μ⁺μ⁻→e⁺e⁻ννννγ via the decay of a muon
given a photon energy `eᵧ`, a center of mass energy `e_cm`, the scalar mediator
model `mod`.
"""
dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator) =
    2 * dndeᵧ_μ_decay(eᵧ, e_cm / 2)

"""
    dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::AbstractScalarMediator)

Compute the photon spectrum from χχ̄ → μ⁺μ⁻→e⁺e⁻ννννγ via the decay of a muon
given a photon energy `eᵧ`, a center of mass energy `e_cm`, the HeavyQuark
model `mod`. Note this is zero since there are no coupling between the scalar
mediator and leptons in the HeavyQuark model.
"""
dndeᵧ_χχ_to_μμ(eᵧ::Real, e_cm::Real, mod::HeavyQuark) = zero(eᵧ)
