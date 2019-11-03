"""
    positron_spectrum_charged_pion(ep::Real, eπ::Real)

Returns the electron/positron spectrum at an electron/positron energy `ep` from
the charged-pion given an arbitrary charged-pion energy `eπ`.
"""
function positron_spectrum_charged_pion(ep::Real, eπ::Real)
    eπ < mπ && return 0.0

    # Line contributions from π → e ν
    line = boosted_delta_function(eπ, mπ, ep, me, (mπ^2 + me^2) / (2 * mπ))
    # Continuum spectrum from π → μν -> e ν ν ν
    spectrum_rf(ep) = positron_spectrum_muon(ep, (mπ^2 + mμ^2) / (2 * mπ))
    continuum = boost_spectrum(spectrum_rf, eπ, mπ, ep, me)

    br_π_eν * line + br_π_μν * continuum
end
