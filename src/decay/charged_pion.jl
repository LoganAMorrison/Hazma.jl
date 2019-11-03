const ENG_GAM_MAX_MU_RF = (mμ^2 - me^2) / (2mμ)
const ENG_MU_PI_RF = (mπ^2 + mμ^2) / (2mπ)

"""
    dnde_pi_to_lnug(eγ::Real, ml::Real)

Return the γ-ray spectrum at photon energy `eγ` from the decay of a charged
pion into a lepton with mass `ml`.
"""
function dnde_pi_to_lnug(eγ::Real, ml::Real)

    x = 2 * eγ / mπ
    r = (ml / mπ)^2

    (x < 0 || (1 - r) < x) && return zero(typeof(eγ))

    # Account for energy-dependence of vector form factor
    F_V = VECTOR_FORM_FACTOR_PI * (1 + VECTOR_FORM_FACTOR_SLOPE_PI * (1 - x))

    # Numerator terms with no log
    f = (r + x - 1) *
        (mπ^2 * x^4 * (AXIAL_FORM_FACTOR_PI^2 + F_V^2) *
         (r^2 - r * x + r - 2 * (x - 1)^2) -
         12 * sqrt(2) * fπ * mπ * r *
         (x - 1) * x^2 * (AXIAL_FORM_FACTOR_PI * (r - 2 * x + 1) + F_V * x) -
         24 * fπ^2 * r * (x - 1) *
         (4 * r * (x - 1) + (x - 2)^2))

    # Numerator terms with log
    g = 12 * sqrt(2) * fπ * r * (x - 1)^2 *
        log(r / (1 - x)) *
        (mπ * x^2 *
         (AXIAL_FORM_FACTOR_PI * (x - 2 * r) - F_V * x) +
         sqrt(2) * fπ *
         (2 * r^2 - 2 * r * x - x^2 + 2 * x - 2))

    αem * (f + g) /
    (24 * π * mπ * fπ^2 * (r - 1)^2 *
     (x - 1)^2 * r * x)
end


"""
    eγ_max(eπ::Real)

Returns the maximum allowed γ-ray energy from a charged pion decay with
energy `eπ` in the laboratory frame.

# Notes
This is computed using the fact that in the mu restframe, the maximum allowed
value is
    eγ_max_μ_rf = (mμ^2 - me^2) / (2mμ)
Then, boosting into the pion rest frame, then to the mu rest frame, we get the
maximum allowed energy in the lab frame.
"""
function eγ_max(eπ::Real)
    βπ = sqrt(1 - (mπ / eπ)^2)
    γπ = eπ / mπ

    βμ = sqrt(1 - (mμ / ENG_MU_PI_RF)^2)
    γμ = ENG_MU_PI_RF / mμ

    ENG_GAM_MAX_MU_RF * γμ * γπ * (1 + βπ) * (1 + βμ)
end

"""
    integrand(cosθ::Real, eγ::Real, eπ::Real[;mode::String="total"])

Returns the integrand of the differential radiative decay spectrum for
the charged pion given a photon angle `cosθ`, a photon energy `eγ`, pion mass
`eπ` given the mode `mode`.
"""
function dnde_pi_integrand(
    cosθ::Real,
    eγ::Real,
    eπ::Real;
    mode::String = "total",
)
    βπ = sqrt(1 - (mπ / eπ)^2)
    γπ = eπ / mπ

    βμ = sqrt(1 - (mμ / ENG_MU_PI_RF)^2)
    γμ = ENG_MU_PI_RF / mμ

    eγ_π_rf = eγ * γπ * (1 - βπ * cosθ)
    jac = 1 / (2γπ * abs(1 - βπ * cosθ))

    dnde_munu = zero(typeof(eγ))
    dnde_munug = zero(typeof(eγ))
    dnde_enug = zero(typeof(eγ))

    if 0 < eγ_π_rf && eγ_π_rf < ENG_GAM_MAX_MU_RF * γμ * (1 + βμ)
        dnde_munu = br_π_μμ * jac *
                    decay_spectrum_muon(eγ_π_rf, ENG_MU_PI_RF)
    end

    dnde_munug = br_π_μμ * jac * dnde_pi_to_lnug(eγ_π_rf, mμ)
    dnde_enug = br_π_eν * jac * dnde_pi_to_lnug(eγ_π_rf, me)

    if mode == "total"
        return dnde_munu + dnde_munug + dnde_enug
    elseif mode == "μν"
        return dnde_munu
    elseif mode == "μνγ"
        return dnde_munug
    elseif mode == "eνγ"
        return dnde_enug
    else
        return zero(typeof(eγ))
    end
end

"""
    decay_spectrum_charged_pion(eγ::Real, eπ::Real[;mode::String=total])

Returns the radiative spectrum value from charged pion given a γ-ray energy
`eγ` and charged pion energy `eπ`. The mode for the decay can be set using
`mode` with options "total", "μν", "μνγ" or "eνγ"
"""
function decay_spectrum_charged_pion(eγ::Real, eπ::Real; mode::String = "total")
    eπ < mπ && return zero(typeof(eγ))
    spectrum_rf(cosθ) = dnde_pi_integrand(cosθ, eγ, eπ; mode = mode);
    quadgk(spectrum_rf, -1, 1)[1]
end
