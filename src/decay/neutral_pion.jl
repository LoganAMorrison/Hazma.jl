"""
    decay_spectrum_neutral_pion(eγ::Real, eπ::Real)

Returns the radiative decay spectrum value from neutral pion given a γ-ray
energy `eγ` and neutral pion energy `eπ`.
"""
function decay_spectrum_neutral_pion(eγ::Real, eπ::Real)
    eπ < NEUTRAL_PION_MASS && return zero(typeof(eγ))
    br_π⁰_γγ * boosted_delta_function(eπ, mπ, eγ, zero(typeof(eγ)), mπ / 2)
end
