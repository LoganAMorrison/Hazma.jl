"""
    dndeₑ_μ_decay_rf(eng_p::Real)

Returns the electron/positron spectrum from a muon in the muon rest frame given
a electron/positron energy `ep`. Note this spectrum is from μ→eνν.
"""
function dndeₑ_μ_decay_μrf(ep::Real)
    r = me / mμ
    ϵ = ep / mμ

    (ϵ < r || (1 + r^2) / 2 < ϵ) && return 0.0

    ((16 * sqrt(ϵ^2 - r^2) * (ϵ * (-3 + 4 * ϵ) + (2 - 3 * ϵ) * r^2)) /
     (mμ * (-1 + 8 * r^2 - 8 * r^6 + r^8 + 24 * r^4 * log(r))))
end

"""
    dndeₑ_μ_decay(ep::Real, eμ::Real)

Returns the electron/positron spectrum at an electron/positron energy `ep` from
the muon given an arbitrary muon energy `eμ`.
"""
function dndeₑ_μ_decay(ep::Real, eμ::Real)
    ed_lb = me
    ed_ub = (mμ^2 + me^2) / (2mμ)
    boost_spectrum(dndeₑ_μ_decay_μrf, eμ, mμ, ep, me; ed_ub=ed_ub, ed_lb=ed_lb)
end
