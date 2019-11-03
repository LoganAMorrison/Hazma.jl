function j_plus(y::Real)
    yconj = 1 - y
    r = (me / mμ)^2
    preFactor = αem * yconj / 6π
    term1 = 3log(yconj / r) - 17 / 2
    term2 = -3log(yconj / r) + 7
    term3 = 2log(yconj / r) - 13 / 3
    preFactor * (term1 + term2 * yconj + term3 * yconj^2)
end

function j_minus(y::Real)
    yconj = 1 - y
    r = (me / mμ)^2
    preFactor = αem * yconj^2 / 6π
    term1 = 3log(yconj / r) - 93 / 12
    term2 = -4log(yconj / r) + 29 / 3
    term3 = 2log(yconj / r) - 55 / 12
    preFactor * (term1 + term2 * yconj + term3 * yconj^2)
end

function dBdy(y::Real)
    (y < 0 || y > 1 - (me / mμ)^2) && return zero(typeof(y))
    (2 / y) * (j_plus(y) + j_minus(y))
end

function decay_spectrum_muon(eγ::Real, eμ::Real)
    spectrum_rf(e) = 2 / mμ * dBdy(2 * e / mμ)
    boost_spectrum(spectrum_rf, eμ, mμ, eγ, zero(typeof(eγ)))
end
