# Generic functions for performing boost integrals.

"""
    boost_energy(ep, mp, ed, md, cosθ)

Returns the boosted energy of the daugther particle when boost from the
lab-frame to the rest-frame of the parent particle.

# Arguments
-`ep::Real`: Energy of parent particle in lab-frame.
-`mp::Real`: Mass of parent particle.
-`ed::Real`: Energy of daughter particle in lab-frame.
-`md::Real`: Mass of daughter particle.
-`cosθ::Real`: Cosine of angle daughter particle makes w.r.t. z-axis in lab-frame
"""
function boost_energy(ep::Real, mp::Real, ed::Real, md::Real, cosθ::Real)
    momentum = sqrt(ep^2 - mp^2)
    γ = ep / mp
    β = momentum / ep

    γ * (ed - β * cosθ * sqrt(ed^2 - md^2))
end

"""
    boost_jacobian(ep, mp, ed, md, cosθ)

Returns the Jacobian for boost integrals when boosting from the lab frame
to the parent particle's rest frame.

# Notes
The Jacobian is given by:
    J = det([∂ER/∂EL ∂ER/∂cosθL; ∂cosθR/∂El ∂cosθR/∂cosθL])
where `ER` is the energy of the daughter particle in the parent particle's
rest-frame, `cosθR` is the cosine of the angle the daughter particle makes
w.r.t. the z-axis. The quantities with `L` are in the lab-frame.

# Arguments
-`ep::Real`: Energy of parent particle in lab-frame.
-`mp::Real`: Mass of parent particle.
-`ed::Real`: Energy of daughter particle in lab-frame.
-`md::Real`: Mass of daughter particle.
-`cosθ::Real`: Cosine of angle daughter particle makes w.r.t. z-axis in lab-frame
"""
function boost_jacobian(ep::Real, mp::Real, ed::Real, md::Real, cosθ::Real)
    momentum = sqrt(ep^2 - mp^2)
    γ = ep / mp
    β = momentum / ep

    ((sqrt(ed^2 - md^2) * (1 - β^2) * γ^2) /
     sqrt(-md^2 + (ed - cosθ * sqrt(ed^2 - md^2) * β)^2 * γ^2))
end

"""
    boosted_delta_function(ep, mp, ed, md, e0)

Boost a δ-function spectrum centered at `e0` from the rest-frame of the
parent particle to the lab-frame.

# Arguments
-`ep::Real`: Energy of parent particle in lab-frame.
-`mp::Real`: Mass of parent particle.
-`ed::Real`: Energy of daughter particle in lab-frame.
-`md::Real`: Mass of daughter particle.
-`e0::Real`: Center of the δ-function spectrum in rest-frame.
"""
function boosted_delta_function(
    ep::Real,
    mp::Real,
    ed::Real,
    md::Real,
    e0::Real,
)
    momentum = sqrt(ep^2 - mp^2)
    γ = ep / mp
    β = momentum / ep

    cond1 = e0 - ed * γ + sqrt(ed^2 - md^2) * β * γ
    cond2 = -e0 + ed * γ + sqrt(ed^2 - md^2) * β * γ

    if cond1 > 0 && cond2 > 0
        return 1 / (2 * sqrt(e0^2 - md^2) * β * γ)
    else
        return 0
    end
end

"""
    boost_spectrum(spectrum_rf, ep, mp, ed, md)

Boost the spectrum of a daughter particle in particle particles rest-frame into the lab-frame.

# Arguments
-`spectrum_rf::Function`: Spectrum in parent particle's rest-frame.
-`ep::Real`: Energy of parent particle in lab-frame.
-`mp::Real`: Mass of parent particle.
-`ed::Real`: Energy of daughter particle in lab-frame.
-`ed_ub::Real`: Upper bound on daughter energy in rest-frame.
-`ed_lw::Real`: Lower bound on daughter energy in rest-frame.
"""
function boost_spectrum(
    spectrum_rf::Function,
    ep::Real,
    mp::Real,
    ed::Real,
    md::Real;
    ed_ub::Real = Inf,
    ed_lb::Real = -Inf,
)
    ep < mp && return zero(typeof(ed))
    ep == mp && return spectrum_rf(ed)
    function integrand(cosθ)
        boostedeng = boost_energy(ep, mp, ed, md, cosθ)
        boostedeng < ed_lb || ed_ub < boostedeng && return zero(typeof(cosθ))
        boost_jacobian(ep, mp, ed, md, cosθ) * spectrum_rf(boostedeng) / 2
    end

    momentum = sqrt(ep^2 - mp^2)
    γ = ep / mp
    β = momentum / ep
    # integration bounds on cosθ
    lb = max((ed - ed_ub / γ) / (β * sqrt(ed^2 - md^2)), -1)
    ub = min((ed - ed_lb / γ) / (β * sqrt(ed^2 - md^2)), 1)
    (lb > 1 || ub < -1) && return zero(typeof(ed))

    quadgk(integrand, lb, ub)[1]
end
