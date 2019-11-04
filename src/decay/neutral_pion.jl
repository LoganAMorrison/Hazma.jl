"""
    dndeᵧ_π⁰_decay(eγ::Real, eπ::Real)

Returns the radiative decay spectrum value from neutral pion given a γ-ray
energy `eγ` and neutral pion energy `eπ`.
"""
function dndeᵧ_π⁰_decay(eγ::Real, eπ::Real)
    eπ < mπ⁰ && return zero(typeof(eγ))
    br_π⁰_γγ * boosted_delta_function(eπ, mπ⁰, eγ, zero(typeof(eγ)), mπ⁰ / 2)
end
