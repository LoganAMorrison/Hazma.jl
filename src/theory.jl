"""A sub-GeV dark matter theory."""
abstract type AbstractTheory end

"""A theory width an s-channel mediator."""
abstract type AbstractTheoryMediator <: AbstractTheory end

# ------------------- #
# Interface functions #
# ------------------- #

"""
    list_annihilation_final_states(mod::AbstractTheory)

Returns a set of strings specifying the available final states for DM
self-annihilation.
"""
list_annihilation_final_states(mod::AbstractTheory) = throw(ErrorException("not implemented"))

"""
    list_decay_final_states(mod::AbstractTheoryMediator)

Returns a set of strings specifying the final states into which the mediator can
decay.
"""
list_decay_final_states(mod::AbstractTheoryMediator) = throw(ErrorException("not implemented"))

"""
    σ_χχ(e_cm, mod, fs)

Returns the cross section for DM annihilating with center-of-mass energy `e_cm`
into the final state `fs`.
"""
σ_χχ(e_cm::Real, mod::AbstractTheory, fs::String) = throw(ErrorException("not implemented"))

"""
    Γ_med(mod, fs)

Returns the mediator partial width for decay into the final state `fs`.
"""
Γ_med(mod::AbstractTheoryMediator, fs::String) = throw(ErrorException("not implemented"))

# TODO: not sure if we'll need this
# σ_χχ_fns(mod::AbstractTheory) = throw(ErrorException("not implemented"))

dndeᵧ_χχ(eᵧ::Real, e_cm::Real, mod::AbstractTheory) = throw(ErrorException("not implemented"))
lines_γ(e_cm::Real, mod::AbstractTheory) = throw(ErrorException("not implemented"))

dndeₑ_χχ(eₑ::Real, e_cm::Real, mod::AbstractTheory) = throw(ErrorException("not implemented"))
lines_e(e_cm::Real, mod::AbstractTheory) = throw(ErrorException("not implemented"))

# --------------------- #
# Implemented functions #
# --------------------- #

"""
    σ_χχ(e_cm, mod)

Returns the total cross section for DM annihilating with center-of-mass energy
`e_cm`.
"""
function σ_χχ(e_cm::Real, mod::AbstractTheory)
    σ_tot = 0.0
    for fs in list_annihilation_final_states(mod)
        σ_tot += σ_χχ(e_cm, mod, fs)
    end
    return σ_tot
end

"""
    br_χχ(e_cm, mod)

Returns a `Dict` of branching fractions for DM annihilating with center-of-mass
energy `e_cm` into each available final state.
"""
function br_χχ(e_cm::Real, mod::AbstractTheory)
    brs = Dict{String,Float64}()

    σ_tot = 0.0
    for fs in list_annihilation_final_states(mod)
        σ_fs = σ_χχ(e_cm, mod, fs)
        brs[fs] = σ_fs
        σ_tot += σ_fs
    end

    if σ_tot != 0
        for (fs, σ_fs) in brs
            brs[fs] = σ_fs / σ_tot
        end
    end

    return brs
end

"""
    Γ_med(mod)

Returns the total mediator decay width.
"""
function Γ_med(mod::AbstractTheoryMediator)
    Γ_tot = 0.0
    for fs in list_decay_final_states(mod)
        Γ_tot += Γ_med(mod, fs)
    end
    return Γ_tot
end

"""
    br_med(e_cm, mod)

Returns a `Dict` of branching fractions for mediator decay into each available
final state.
"""
function br_med(mod::AbstractTheoryMediator)
    brs = Dict{String,Float64}()

    Γ_tot = 0.0
    for fs in list_annihilation_final_states(mod)
        Γ_fs = Γ_χχ(mod, fs)
        brs[fs] = Γ_fs
        Γ_tot += Γ_fs
    end

    if Γ_tot != 0
        for (fs, Γ_fs) in brs
            brs[fs] = Γ_fs / Γ_tot
        end
    end

    return brs
end
