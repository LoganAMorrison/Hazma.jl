"""
    dndeγ_π_lνγ_πrf(eγ::Real, ml::Real)

Return the γ-ray spectrum at photon energy `eγ` from the decay of a
charged-pion into a lepton with mass `ml` in the charged-pion's rest-frame.
"""
function dndeγ_π_lνγ_πrf(eγ::Real, ml::Real)

    x = 2 * eγ / mπ
    r = (ml / mπ)^2

    (x < 0 || (1 - r) < x) && return zero(typeof(eγ))

    # Account for energy-dependence of vector form factor
    F_V = VECTOR_FORM_FACTOR_PI * (1 + VECTOR_FORM_FACTOR_SLOPE_PI * (1 - x))

    # Numerator terms with no log
    f = (r + x - 1) * (mπ^2 * x^4 * (AXIAL_FORM_FACTOR_PI^2 + F_V^2) *
         (r^2 - r * x + r - 2 * (x - 1)^2) -
         12 * sqrt(2) * fπ * mπ * r * (x - 1) * x^2 *
         (AXIAL_FORM_FACTOR_PI * (r - 2 * x + 1) + F_V * x) -
         24 * fπ^2 * r * (x - 1) * (4 * r * (x - 1) + (x - 2)^2))

    # Numerator terms with log
    g = 12 * sqrt(2) * fπ * r * (x - 1)^2 * log(r / (1 - x)) *
        (mπ * x^2 * (AXIAL_FORM_FACTOR_PI * (x - 2 * r) - F_V * x) +
         sqrt(2) * fπ * (2 * r^2 - 2 * r * x - x^2 + 2 * x - 2))

    αem * (f + g) / (24 * π * mπ * fπ^2 * (r - 1)^2 * (x - 1)^2 * r * x)
end

"""
    dndeγ_π_decay(eγ::Real, eπ::Real[;mode::String=total])

Returns the radiative decay spectrum from charged pion given a γ-ray energy
`eγ` and charged pion energy `eπ`. The mode for the decay can be set using
`mode` with options "total", "μν", "μνγ" or "eνγ"
"""
function dndeγ_π_decay(eγ::Real, eπ::Real; mode::String = "total")
    eπ < mπ && return zero(typeof(eγ))

    dndeγ_π_eνγ_πrf(engγ) = br_π_eν * dndeγ_π_lνγ_πrf(engγ, me)
    dndeγ_π_μνγ_πrf(engγ) = br_π_μν * dndeγ_π_lνγ_πrf(engγ, mμ)
    dndeγ_π_μν_πrf(engγ) =
        br_π_μν * dndeγ_μ_decay(engγ, (mπ^2 + mμ^2) / (2 * mπ))
    spectrum_rf(engγ) =
        (dndeγ_π_eνγ_πrf(engγ) + dndeγ_π_μνγ_πrf(engγ) + dndeγ_π_μν_πrf(engγ))

    spec = zero(typeof(eγ))
    if mode == "total"
        return boost_spectrum(spectrum_rf, eπ, mπ, eγ, 0.0)
    elseif mode == "μν"
        return boost_spectrum(dndeγ_π_μν_πrf, eπ, mπ, eγ, 0.0)
    elseif mode == "μνγ"
        return boost_spectrum(dndeγ_π_μνγ_πrf, eπ, mπ, eγ, 0.0)
    elseif mode == "eνγ"
        return boost_spectrum(dndeγ_π_eνγ_πrf, eπ, mπ, eγ, 0.0)
    else
        return zero(typeof(eγ))
    end
end
