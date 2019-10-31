function j_plus(y::Real)
    yconj = 1 - y
    r = (ELECTRON_MASS / MUON_MASS)^2
    preFactor = ALPHA_EM * yconj / 6π
    term1 = 3log(yconj / r) - 17 / 2
    term2 = -3log(yconj / r) + 7
    term3 = 2log(yconj / r) - 13 / 3
    preFactor * (term1 + term2 * yconj + term3 * yconj^2)
end

function j_minus(y::Real)
    yconj = 1 - y
    r = (ELECTRON_MASS / MUON_MASS)^2
    preFactor = ALPHA_EM * yconj^2 / 6π
    term1 = 3log(yconj / r) - 93 / 12
    term2 = -4log(yconj / r) + 29 / 3
    term3 = 2log(yconj / r) - 55 / 12
    preFactor * (term1 + term2 * yconj + term3 * yconj^2)
end

function dBdy(y::Real)
    (y < 0 || y > 1 - (ELECTRON_MASS / MUON_MASS)^2) && return zero(typeof(y))
    (2.0 / y) * (j_plus(y) + j_minus(y))
end

function dnde_muon_integrand(cl::Real, eγ::Real, eμ::Real)
    β = sqrt(1 - (MUON_MASS / eμ)^2)
    γ = MUON_MASS / eμ
    eγ_μrf = γ * eγ * (1 - β * cl)
    dBdy((2 / mμ) * eγ_μrf) / (eμ * (1 - cl * β))
end

function decay_spectrum_muon(eγ::Real, eμ::Real)
    eμ < mμ && return zero(typeof(eγ))

    β = sqrt(1 - (MUON_MASS / eμ)^2)
    γ = MUON_MASS / eμ

    eγ_max = (mμ - me^2 / mμ) * γ * (1 + β) / 2

    (eγ < 0 || eγ > eγ_max) && return zero(typeof(eγ))

    nodes, weights = gausslegendre(15)
    sum(weights .* dnde_muon_integrand.(nodes, eγ, eμ))
end
