
"""
    positron_spectrum_muon_rf(eng_p::Real)

Returns the electron/positron spectrum from a muon in the muon rest frame given
a electron/positron energy `ep`.
"""
function positron_spectrum_muon_rf(ep::Real)
    r = me / mμ
    s = me^2 - 2 * ep * mμ + mμ^2
    smax = (mμ - me)^2
    smin = 0.0

    (s <= smin || smax <= s) && return 0.0

    return 2 * mμ *
           (2 * (mμ^4 * (r^2 - 1)^2 + mμ^2 * (1 + r^2) * s - 2 * s^2) *
            sqrt(mμ^4 * (r^2 - 1)^2 - 2 * mμ^2 * (1 + r^2) * s + s^2)) / mμ^8
end

"""
    positron_spectrum_muon_integrand(cosθ::Real, ep::Real, eμ::Real)

Returns the integrand for the boost integral used to compute the electron/
positron spectrum from the muon given a electron/positron angle `cosθ`, energy
`ep` and muon energy `eμ`.
"""
function positron_spectrum_muon_integrand(cosθ::Real, ep::Real, eμ::Real)
    ep < me && return 0.0

    p = sqrt(ep^2 - me^2)
    γ = eμ / mμ
    β = sqrt(1.0 - (mμ / eμ)^2)
    eμrf = γ * (ep - p * β * cosθ)
    jac = (p / (2.0 *
            sqrt((1 + (β * cosθ)^2) * ep^2 - (1 + β^2 * (-1 + cosθ^2)) * me^2 -
                 2 * β * cosθ * ep * p) *
            γ))

    positron_spectrum_muon_rf(eμrf) * jac
end

"""
    positron_spectrum_muon(ep::Real, eμ::Real)

Returns the electron/positron spectrum at an electron/positron energy `ep` from
the muon given an arbitrary muon energy `eμ`.
"""
function positron_spectrum_muon(ep::Real, eμ::Real)
    eμ < mμ && return 0.0

    nodes, weights = gausslegendre(15)
    sum(weights .* positron_spectrum_muon_integrand.(nodes, ep, eμ))
end
