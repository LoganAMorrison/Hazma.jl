# functions for computing the γ-ray spectrum from  muon decays: μ → e ν ν γ

"""
    j_plus(y::Real)

Form factor in differential branching fraction of radiative muon decay given
a scaled photon energy `y` = 2Eᵧ / mμ. See p.18, eqn (54) of
arXiv:hep-ph/9909265.
"""
function j_plus(y::Real)
    yconj = 1 - y
    r = (me / mμ)^2
    preFactor = αem * yconj / 6π
    term1 = 3 * log(yconj / r) - 17 / 2
    term2 = -3 * log(yconj / r) + 7
    term3 = 2 * log(yconj / r) - 13 / 3
    preFactor * (term1 + term2 * yconj + term3 * yconj^2)
end

"""
    j_minus(y::Real)

Form factor in differential branching fraction of radiative muon decay given
a scaled photon energy `y` = 2Eᵧ / mμ. See p.18, eqn (55) of
arXiv:hep-ph/9909265.
"""
function j_minus(y::Real)
    yconj = 1 - y
    r = (me / mμ)^2
    preFactor = αem * yconj^2 / 6π
    term1 = 3 * log(yconj / r) - 93 / 12
    term2 = -4 * log(yconj / r) + 29 / 3
    term3 = 2 * log(yconj / r) - 55 / 12
    preFactor * (term1 + term2 * yconj + term3 * yconj^2)
end

"""
    decay_spectrum_muon_rf(eγ)

Returns the radiative decay spectrum from muon given a γ-ray energy
`eγ` in the muon rest-frame.
"""
function decay_spectrum_muon_rf(eγ::Real)
    y = 2 * eγ / mμ
    (y < 0 || y > 1 - (me / mμ)^2) && return zero(typeof(y))
    2 / mμ * (2 / y) * (j_plus(y) + j_minus(y))
end

"""
    decay_spectrum_muon(eγ, eμ)

Returns the radiative decay spectrum from muon given a γ-ray energy
`eγ` and muon energy `eμ`.
"""
decay_spectrum_muon(eγ::Real, eμ::Real) =
    boost_spectrum(decay_spectrum_muon_rf, eμ, mμ, eγ, zero(typeof(eγ)))
