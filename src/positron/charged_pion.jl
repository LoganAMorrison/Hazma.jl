"""
    dndeₑ_π_decay(ep::Real, eπ::Real)

Returns the electron/positron spectrum at an electron/positron energy `ep` from
the charged-pion given an arbitrary charged-pion energy `eπ`.
"""
function dndeₑ_π_decay(ep::Real, eπ::Real)
    eπ < mπ && return 0.0

    # compute bounds on integration
    eμ_πrf = (mπ^2 + mμ^2) / (2 * mπ)
    γμ = eμ_πrf / mμ
    βμ = sqrt(1 - mμ^2 / eμ_πrf^2)
    ub = γμ * (me^2 + mμ^2 + βμ * sqrt((me^2 + mμ^2)^2 - 4 * me^2 * mμ^2)) /
         (2 * mμ)

    # Line contributions from π → e ν
    line = boosted_delta_function(eπ, mπ, ep, me, (mπ^2 + me^2) / (2 * mπ))
    # Continuum spectrum from π → μν -> e ν ν ν
    spectrum_rf(ep) = dndeₑ_μ_decay(ep, (mπ^2 + mμ^2) / (2 * mπ))
    continuum = boost_spectrum(
        spectrum_rf,
        eπ,
        mπ,
        ep,
        me;
        ed_ub = ub,
        ed_lb = me,
    )

    br_π_eν * line + br_π_μν * continuum
end
