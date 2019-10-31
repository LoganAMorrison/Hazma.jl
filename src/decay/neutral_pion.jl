"""
    decay_spectrum_neutral_pion(eγ::Real, eπ::Real)

Returns the radiative decay spectrum value from neutral pion given a γ-ray
energy `eγ` and neutral pion energy `eπ`.
"""
function decay_spectrum_neutral_pion(eγ::Real, eπ::Real)
    eπ < NEUTRAL_PION_MASS && return zero(typeof(eγ))

    β = sqrt(1 - (NEUTRAL_PION_MASS / eπ)^2)
    if eπ * (1 - β) / 2 <= eγ <= eπ * (1 + β) / 2
        return BR_PI0_TO_GG * 2 / (eπ * β)
    else
        return zero(typeof(eγ))
    end
end
