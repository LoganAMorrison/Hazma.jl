"""
    positron_spectrum_charged_pion_integrand(cosθ::Real, ep::Real, eπ::Real)

Returns the integrand for the boost integral used to compute the electron/
positron spectrum from the charged pion given a electron/positron angle
`cosθ`, energy `ep` and muon energy `eμ`.
"""
function positron_spectrum_charged_pion_integrand(
    cosθ::Real, ep::Real, eπ::Real)

    ep < me && return 0.0

    p = sqrt(ep^2 - me^2)
    γ = eπ / mπ
    β = sqrt(1.0 - (mπ / eπ)^2)
    jac = (p / (2. * sqrt((1 + (β * cosθ)^2) * ep^2 -
           (1 + β^2 * (-1 + cosθ^2)) * me^2 -
            2 * β * cosθ * ep * p) * γ))

    epπrf = γ * (ep - p * β * cosθ)
    eμπrf = (mπ^2 + mμ^2) / (2 * mπ)

    BR_PI_TO_MUNU * jac * positron_spectrum_muon(epπrf, eμπrf)
end

"""
    positron_spectrum_charged_pion(ep::Real, eπ::Real)

Returns the electron/positron spectrum at an electron/positron energy `ep` from
the charged-pion given an arbitrary charged-pion energy `eπ`.
"""
function positron_spectrum_charged_pion(ep::Real, eπ::Real)
    eπ < mπ && return 0.0

    nodes, weights = gausslegendre(15)
    sum(weights .* positron_spectrum_charged_pion_integrand.(nodes, ep, eπ))
end
