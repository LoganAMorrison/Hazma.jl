"""
    positive(x::Real)

Return `true` if `x` is greater than zero.
"""
positive(x::Real) = x > 0

"""
    nonnegative(x::Real)

Return `true` if `x` is greater than or equal to zero.
"""
nonnegative(x) = x ≥ 0


"""
  horner(x, p...)

Convert an expression of the form a₀ + a₁x + a₂x² + a₃x³ + ... into Horner form:
  a₀ + a₁x + a₂x² -> a₀ + x(a₁ + x(a₂ + x(a₃ + ⋅⋅⋅)
and combines all expressions of the form `ax + y` operatations into `muladd`.

# Example
horner(x, 1, 2, 3) = 1 + 2x + 3x² -> muladd(x, muladd(3,x,2), 1)
"""
macro horner(x, p...)
    ex = esc(p[end])
    for i = length(p)-1:-1:1
        ex = :(muladd(t, $ex, $(esc(p[i]))))
    end
    Expr(:block, :(t = $(esc(x))), ex)
end

# ------------------------------------------------------------------- #
# Pure Julia implementations of modified bessel functions of 2nd kind #
# ------------------------------------------------------------------- #

"""
	besselk0_small_x_pade_logx(z::Real)

Pade approximant of small `z = x^2` approximation with coefficient log(x) for
the modified bessel function of order 0.
"""
function besselk0_small_x_pade_logx(z::Real)
    num = @horner(
        z,
        -1.0,
        -0.24138285628218,
        -0.0135080396880856,
        -0.000308608673200753,
        -3.59813777913474e-6,
        -2.39537920873459e-8,
        -9.64782224098546e-11,
        -2.39973630994234e-13,
        -3.6220564063639e-16,
        -3.06341714046988e-19,
        -1.12488231252247e-22,
    )
    den = @horner(
        z,
        1.0,
        -0.00861714371782045,
        0.000037325617540672,
        -1.07638371248223e-7,
        2.30308733362263e-10,
        -3.84742171501404e-13,
        5.12223208206566e-16,
        -5.41420521800473e-19,
        4.39231123172861e-22,
        -2.4982196441535e-25,
        7.66670699554804e-29,
    )
    num / den
end

"""
	besselk0_small_x_pade(z::Real)

Pade approximant of small x approximation without logs for the modified bessel
function of order 0.
"""
function besselk0_small_x_pade(z::Real)
    num = @horner(
        z,
        0.11593151565841242,
        0.2779669968140355,
        0.02280874457582517,
        0.0006355510644535478,
        8.445064385629965e-6,
        6.190852691234865e-8,
        2.689915420106763e-10,
        7.120509070871071e-13,
        1.1328286076478593e-15,
        1.002684743754697e-18,
        3.831825269223243e-22,
    )
    den = @horner(
        z,
        1.0,
        -0.008762777703698878,
        0.000038640008315635896,
        -1.1357363845585156e-7,
        2.4802453234587166e-10,
        -4.2354464596974406e-13,
        5.774298134837363e-16,
        -6.262762197904977e-19,
        5.225648318588274e-22,
        -3.065499525380141e-25,
        9.735295128548586e-29,
    )
    num / den
end

"""
	besselk0_large_x_pade(z)

Pade approximant of large x approximation with e^-x/√x factored out for the
modified bessel function of the 2nd kind of order 0.
"""
function besselk0_large_x_pade(x::Real)
    num = @horner(
        x,
        1.2533141373155,
        141.331887132623,
        6942.55823845144,
        196059.379528513,
        3.53857474769014e6,
        4.29329824457297e7,
        3.58713372871623e8,
        2.08024894901955e9,
        8.33492502402265e9,
        2.26975694930578e10,
        4.076523942544e10,
        4.59887295747372e10,
        3.0123468462625e10,
        1.00104935562229e10,
        1.29061584510892e9,
        2.97679409536703e7,
    )
    den = @horner(
        x,
        1.0,
        112.891530692253,
        5553.40115977167,
        157119.062627797,
        2.84263172400885e6,
        3.46002396834131e7,
        2.90347929369435e8,
        1.69385095851949e9,
        6.84387235591249e9,
        1.88643752101451e10,
        3.45007614511026e10,
        4.00413944620127e10,
        2.74927582795927e10,
        9.95043548034518e9,
        1.53544999786975e9,
        6.12174527802927e7,
    )
    num / den
end

"""
	besselk1_small_x_pade_logx(z::Real)

Pade approximant of small `z = x^2` approximation with coefficient log(x) for
the modified bessel function of 2nd kind of order 1.
"""
function besselk1_small_x_pade_logx(z::Real)
    num = @horner(
        z,
        0.5,
        0.0584439151701907,
        0.00211364524554443,
        0.0000351446989590961,
        3.18464978980487e-7,
        1.71846370827642e-9,
        5.77767674391542e-12,
        1.2261731998997e-14,
        1.60620368680067e-17,
        1.19522891185578e-20,
        3.90527439361035e-24,
    )
    den = @horner(
        z,
        1.0,
        -0.00811216965961851,
        0.0000329783652078475,
        -8.8958533386728e-8,
        1.77392868020876e-10,
        -2.75069415038282e-13,
        3.38397444674368e-16,
        -3.28865819658875e-19,
        2.43917080106789e-22,
        -1.26027084469828e-25,
        3.48769663012292e-29,
    )
    num / den
end

"""
	besselk0_small_x_pade(x::Real)

Pade approximant of small x approximation without logs for the modified bessel
function of the 2nd kind of order 1.
"""
function besselk1_small_x_pade(x::Real)
    num = @horner(
        x,
        -0.3079657578292062,
        -0.0828319564150824,
        -0.003948914845123864,
        -0.00007714903403507421,
        -7.818995663494835e-7,
        -4.593787841488269e-9,
        -1.6535766309440503e-11,
        -3.714375432350668e-14,
        -5.1073001491694574e-17,
        -3.964303282935162e-20,
        -1.3444699412781958e-23,
    )
    den = @horner(
        x,
        1.0,
        -0.008243654526606038,
        0.00003409079911031547,
        -9.365135738273977e-8,
        1.9042924092208716e-10,
        -3.0153254895785485e-13,
        3.7942153081309085e-16,
        -3.778584837430635e-19,
        2.8781246406305737e-22,
        -1.5310603544084845e-25,
        4.375644955511107e-29,
    )
    num / den
end

"""
	besselk1_large_x_pade(x::Real)

Pade approximant of large x approximation with e^-x/√x factored out for the
modified bessel function of the 2nd kind of order 1.
"""
function besselk1_large_x_pade(x::Real)
    num = @horner(
        x,
        1.2533141373155,
        140.214272883603,
        6831.20156914534,
        191298.484366315,
        3.42376570691759e6,
        4.12068812662786e7,
        3.41844499627513e8,
        1.97205766074477e9,
        7.88811226926147e9,
        2.15834949563827e10,
        3.94055246985218e10,
        4.61577651815379e10,
        3.26603924021509e10,
        1.2671639926028e10,
        2.25882216818188e9,
        1.20800792278632e8,
    )
    den = @horner(
        x,
        1.0,
        111.499803538027,
        5408.81502492662,
        150618.765390688,
        2.67591031695531e6,
        3.18919798142894e7,
        2.61091844605052e8,
        1.47904816548709e9,
        5.76683394861962e9,
        1.52090762843596e10,
        2.62921371781387e10,
        2.83179496473139e10,
        1.75105653887285e10,
        5.40116812236514e9,
        6.26500374809249e8,
        1.15380997710812e7,
    )
    num / den
end


"""
	besselk0(x::Real)

Compute the modified Bessel function of the second kind of order 0 at `x`.
"""
function besselk0(x::Real)
    if x <= 3.0
        z = x^2
        term = besselk0_small_x_pade_logx(z) * log(x)
        return besselk0_small_x_pade(z) + term
    else
        z = 1 / x
        return exp(-x) * besselk0_large_x_pade(z) / sqrt(x)
    end
end

"""
	besselk1(x::Real)

Compute the modified Bessel function of the second kind of order 1 at `x`.
"""
function besselk1(x::Real)
    if x <= 3.0
        z = x^2
        term = besselk1_small_x_pade_logx(z) * log(x)
        return x * (besselk1_small_x_pade(z) + term) + 1 / x
    else
        z = 1 / x
        return exp(-x) * besselk1_large_x_pade(z) / sqrt(x)
    end
end

"""
	besselkn(n::Int, x::Real)

Compute the modified Bessel function of the second kind of order `n` at `x`.
"""
function besselkn(n::Int, x::Real)
    n == 0 && return besselk0(x)
    n == 1 && return besselk1(x)

    tox = 2 / x
    bkm = besselk0(x)
    bk = besselk1(x)
    for j = 1:n-1
        bkp = bkm + j * tox * bk
        bkm = bk
        bk = bkp
    end
    bk
end

"""
	besselk2(x::Real)

Compute the modified Bessel function of the second kind of order 2 at `x`.
"""
besselk2(x::Real) = besselk0(x) + 2 * besselk1(x) / x

"""
	besselk3(x::Real)

Compute the modified Bessel function of the second kind of order 3 at `x`.
"""
besselk3(x::Real) = ((x^2 + 8) * besselk1(x) + 4 * x * besselk0(x)) / x^2
