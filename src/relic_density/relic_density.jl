# Load SM data: T, √g⋆, heff, geff
_sm_data = readdlm(joinpath(string(@__DIR__), "smdof.csv"), ',', skipstart = 1);
_sm_tempetatures = _sm_data[:, 1] * 1e3;  # convert to MeV
_sm_sqrt_gstars = _sm_data[:, 2];
_sm_heff = _sm_data[:, 3];

# Interpolating function for SM's sqrt(g_star)
_sm_sqrt_gstar = LinearInterpolation(
    _sm_tempetatures,
    _sm_sqrt_gstars,
    extrapolation_bc = Flat(),
);
# Interpolating function for SM d.o.f. stored in entropy: h_eff
_sm_heff = LinearInterpolation(
    _sm_tempetatures,
    _sm_heff,
    extrapolation_bc = Flat(),
)


"""
    sm_dof_entropy(T::Real)

Compute the d.o.f. stored in entropy of the Standard Model at a temperature `T`.
"""
sm_dof_entropy(T::Real) = _sm_heff(T);


"""
    sm_sqrt_gstar(T::Real)

Compute the square-root of g-star of the Standard Model at a temperature `T`.
"""
sm_sqrt_gstar(T::Real) = _sm_sqrt_gstar(T);


"""
    sm_entropy_density(T::Real)

Compute the entropy density of the Standard Model at a temperature `T`.
"""
sm_entropy_density(T::Real) = 2 * π^2 / 45 * _sm_heff(T) * T^3;



"""
    neq(T::Real, m::Float64[; g=2.0, stats=:fermion])

Compute the equilibrium number density of a particle with mass
`m` and temperature `T`. For fermions, use `stats=:fermion`
(`stats=:boson` for bosons).
"""
function neq(T::Real, m::Float64; g::Float64 = 2.0, stats::Symbol = :fermion)
    if m == 0
        # if particle is massless, use analytic expression.
        # fermion: 7 / 8 zeta(3) / pi^2
        # boson: zeta(3) / pi^2
        nbar = (stats == :fermion) ? 0.0913453711751798 : 0.121793828233573
    else
        # use sum-over-bessel function representation of neq
        # nbar = x^2 sum_n (\pm 1)^{n+1}/n k_2(nx)
        eta = (stats == :fermion) ? -1 : 1
        x = m / T

        nbar = x^2 * (eta^(1 + 1) / 1 * besselk2(1 * x) +
                eta^(2 + 1) / 2 * besselk2(2 * x) +
                eta^(3 + 1) / 3 * besselk2(3 * x) +
                eta^(4 + 1) / 4 * besselk2(4 * x) +
                eta^(5 + 1) / 5 * besselk2(5 * x)) / (2 * π^2)
    end
    g * nbar * T^3
end


"""
    yeq(T::Real, m::Float64[;g=2.0, stats=:fermion])

Compute the equilibrium value of Y, the comoving number density `neq / s` where
`s` is the SM entropy density.
"""
function yeq(T::Real, m::Float64; g::Float64 = 2.0, stats::Symbol = :fermion)
    neq(T, m; g = g, stats = stats) / sm_entropy_density(T)
end


"""
    weq(T::Real, m::Float64, g::Float64=2.0, stats::Symbol = :fermion)

Compute the equilibrium value of W, the natural log of the comoving number
density `Y` = `neq / s` where `s` is the SM entropy density.
"""
function weq(T::Real, m::Float64; g::Float64 = 2.0, stats::Symbol = :fermion)
    log(neq(T, m; g = g, stats = stats) / sm_entropy_density(T))
end

"""
    thermal_cross_section(x::Real, model::AbstractTheory)

Compute the thermally average cross section for the dark
matter particle of the given model.
"""
function thermal_cross_section(x::T, model::AbstractTheory) where {T<:Real}
    x > 300 && return zero(T)

    pf::T = x / (2 * besselk2(x))^2
    integrand(z::T) =
        (z^2 * (z^2 - 4) * besselk1(x * z) * σ_χχ(model.mχ * z, model))
    pf * quadgk(integrand, 2.0, max(50 / x, 150))[1]
end

# ----------------------------------------------- #
# Functions for solving the Bolzmann equation     #
# see arXiv:1204.3622v3 for detailed explaination #
# ------------------------------- --------------- #

"""
    boltzmann_eqn!(dw, w, logx, model::AbstractTheory)

Compute the RHS of the Boltzmann equation. Here the RHS is given by dW/dlogx,
with W = log(neq / sm_entropy_density). Note that this function is intended to
be used with `DifferentialEquations.jl`.
"""
function boltzmann_eqn!(dw, w, model::AbstractTheory, logx)
    mx = model.mχ
    x = exp(logx)
    T = mx / x
    pf = -sqrt(π / 45) * PLANK_MASS * mx * sm_sqrt_gstar(T) / x
    _weq = weq(T, mx, g = 2.0)
    sv = thermal_cross_section(x, model)

    dw[1] = pf * sv * (exp(w[1]) - exp(2.0 * _weq - w[1]))
end

"""
    jacobian_boltzmann_eqn!(logx, w, model::AbstractTheory)

Compute the Jacobian of the RHS of the Boltzmann equation with
respect to the log of the comoving equilibrium number density.
"""
function jacobian_boltzmann_eqn!(J, w, model::AbstractTheory, logx)
    mx = model.mχ
    x = exp(logx)
    T = mx / x
    pf = -sqrt(π / 45) * PLANK_MASS * mx * sm_sqrt_gstar(T) / x
    _weq = weq(T, mx, g = 2.0)
    sv = thermal_cross_section(x, model)

    J[1, 1] = pf * sv * (exp(w[1]) + exp(2.0 * _weq - w[1]))
end

"""
    solve_boltzmann(model::AbstractTheory)

Solve the Boltzmann equation for the log of the dark matter comoving number
density as a function of `logx` - which is the log of the dark matter mass over
its temperature.
"""
function solve_boltzmann(
    model::AbstractTheory;
    x0 = 1.0,
    xf = nothing,
    alg = radau(),
    reltol = 1e-5,
    abstol = 1e-3,
)
    mx = model.mχ
    T0 = mx / x0
    w0 = [weq(T0, mx; g = 2.0)]
    logx0 = log(x0)
    logxf = (xf == nothing) ? logx0 + 7.0 : log(xf)

    fj = ODEFunction(boltzmann_eqn!; jac = jacobian_boltzmann_eqn!)
    prob = ODEProblem(fj, w0, (logx0, logxf), model)
    solve(prob, alg, reltol = reltol, abstol = abstol)
end


# ----------------------------------------------------------- #
# Functions for computing the relic density semi-analytically #
# see arXiv:1204.3622v3 for detailed explaination             #
# ----------------------------------------------------------- #

"""
    xstar_root_eqn(xstar, model::AbstractTheory; δ)

Returns residual of root equation used to solve for x_star given the current
value `xstar` and the DM model `model` assuming that the DM freezes out when
`Y ~ (1+δ) Yeq` with `δ` by default set to (√5 - 1) / 2.

See Eqn.(14) of arXiv:1204.3622v3 for similar expressions. Note that our
result is more exact since we do not assume `xstar` is large enough that
Yeq ~ x^{3/2} e^{-x} / h_sm. This may cause a bit of a slow down.
"""
function xstar_root_eqn(
    xstar::Real,
    model::AbstractTheory;
    δ::Float64 = 0.5 * (sqrt(5) - 1),
)
    Δ = δ * (2 + δ) / (1 + δ)
    T = model.mχ / xstar
    λ = sqrt(π / 45) * model.mχ * PLANK_MASS * sm_sqrt_gstar(T)
    tcs = thermal_cross_section(xstar, model)
    _yeq = yeq(T, model.mχ)
    dyeq = ForwardDiff.derivative(x -> yeq(model.mχ / x, model.mχ), xstar)
    return xstar^2 * dyeq + λ * Δ * tcs * _yeq^2
end

"""
    compute_xstar(model::AbstractTheory; δ::Float64=0.5 * (sqrt(5) - 1))

Computes to value of `xstar`: the value of dm_mass / temperature such that
the DM begins to freeze out.
"""
function compute_xstar(model::AbstractTheory; δ::Float64 = 0.5 * (sqrt(5) - 1))
    f(x) = xstar_root_eqn(x, model; δ = δ)
    find_zero(f, (0.01, 100.0))
end

"""
    compute_alpha(model, xstar)

Computes the value of the integral of RHS of the Boltzmann equation with
Yeq set to zero from x_{\\star} to x_{\\mathrm{f.o.}}.
"""
function compute_alpha(model::AbstractTheory, xstar)
    pf = sqrt(π / 45) * model.mχ * PLANK_MASS

    integrand(x) =
        (sm_sqrt_gstar(model.mχ / x) * thermal_cross_section(x, model) / x^2)
    pf * quadgk(integrand, xstar, 100 * xstar)[1]
end

"""
    relic_density(model::AbstractTheory; semi_analytic=true, kwargs...)

Solves the Boltzmann equation and returns the relic density
computed from the final dark matter comoving number density.

# Arguments
- `model::AbstractTheory`: Dark matter model.
- `semi_analytic::Bool` : If `True`, the relic density is computed using semi-analytical methods, otherwise the Boltzmann equation is numerically solved.

# kwargs
If `semi_analtyical` is `True`, the accepted kwargs are:
    delta: float, optional
            Value of `delta` assumed for when DM begins to freeze out.
            Default value is the solution to
            delta * (2 + delta) / (1 + delta) = 1, i.e.,
            delta = (sqrt(5) - 1) / 2 = 0.618033988749895. See Eqn.(13) of
            arXiv:1204.3622v3 for details and other used values of delta.
            Value of xstar is logarithmically sensitive to this number.

If `semi_analytical` is `False`, accepted kwargs are:
    x0: float, optional
        Initial value of x = mass / temperature. Default is `1.0`.
    xf: float, optional
        Final value of x = mass / temperature. Default is `1000`
        times the initial starting value.
    method: string, optional
        Method used to solve the Boltzmann equation. Default is
        'Radau'.
    rtol: float, optional
        Relative tolerance used to solve the Boltzmann equation.
        Default is `1e-3`.
    atol: float, optional
        Absolute tolerance used to solve the Boltzmann equation.
        Default is `1e-6`.
"""
function relic_density(model::AbstractTheory; semi_analytic = true, kwargs...)
    if semi_analytic
        δ = (:δ in keys(kwargs)) ? kwargs[:δ] : 0.5 * (sqrt(5) - 1)
        xstar = compute_xstar(model; δ = δ)
        alpha = compute_alpha(model, xstar)
        ystar = yeq(model.mχ / xstar, model.mχ)
        Y0 = ystar / (1 + ystar * alpha)
    else
        x0 = (:x0 in keys(kwargs)) ? kwargs[:x0] : 1.0
        xf = (:xf in keys(kwargs)) ? kwargs[:xf] : nothing
        alg = (:alg in keys(kwargs)) ? kwargs[:alg] : radau()
        reltol = (:reltol in keys(kwargs)) ? kwargs[:reltol] : 1e-5
        abstol = (:abstol in keys(kwargs)) ? kwargs[:abstol] : 1e-3
        sol = solve_boltzmann(
            model;
            x0 = x0,
            xf = xf,
            alg = alg,
            reltol = reltol,
            abstol = abstol,
        )
        Y0 = exp(sol[end][1])
    end
    Y0 * model.mχ * SM_ENTROPY_DENSITY_TODAY / CRITICAL_ENERGY_DENSITY
end
