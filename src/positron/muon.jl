"""
    positron_spectrum_muon_rf(eng_p::Real)

Returns the electron/positron spectrum from a muon in the muon rest frame given
a electron/positron energy `ep`. Note this spectrum is from μ→eνν.
"""
function positron_spectrum_muon_rf(ep::Real)
    r = me / mμ
    ϵ = ep / mμ

    (ϵ < r || (1 + r^2) / 2 < ϵ) && return 0.0

    ((16 * sqrt(ϵ^2 - r^2) * (ϵ * (-3 + 4 * ϵ) + (2 - 3 * ϵ) * r^2)) /
     (mμ * (-1 + 8 * r^2 - 8 * r^6 + r^8 + 24 * r^4 * log(r))))
end

"""
    positron_spectrum_muon_integrand(cosθ::Real, ep::Real, eμ::Real)

Returns the integrand for the boost integral used to compute the electron/
positron spectrum from the muon given a electron/positron angle `cosθ`, energy
`ep` and muon energy `eμ`.
"""
function positron_spectrum_muon_integrand(cosθ::Real, ep::Real, eμ::Real)
    ep < me && return 0.0

    p = sqrt(eμ^2 - mμ^2)
    γ = eμ / mμ
    β = p / eμ
    ϵ = ep / mμ

    eμrf = γ * (ep - sqrt(ep^2 - me^2) * β * cosθ)

    r = me / mμ
    ϵrf = eμrf / mμ

    (1 + r^2) / 2 < ϵrf && return 0.0

    jac = ((((1 - β^2) * γ^2 * sqrt(ϵ^2 - r^2)) /
            sqrt(-r^2 + γ^2 * (ϵ - β * cosθ * sqrt(ϵ^2 - r^2))^2)))

    positron_spectrum_muon_rf(eμrf) * jac / 2
end

"""
    positron_spectrum_muon(ep::Real, eμ::Real)

Returns the electron/positron spectrum at an electron/positron energy `ep` from
the muon given an arbitrary muon energy `eμ`.
"""
function positron_spectrum_muon(ep::Real, eμ::Real)
    eμ < mμ && return 0.0
    eμ == mμ && return positron_spectrum_muon_rf(ep)
    f(cosθ) = positron_spectrum_muon_integrand(cosθ, ep, eμ)
    quadgk(f, -1, 1)[1]
end
