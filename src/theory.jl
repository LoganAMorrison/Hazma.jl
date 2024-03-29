"""A sub-GeV dark matter theory."""
abstract type AbstractTheory end

"""A theory with an s-channel mediator."""
abstract type AbstractTheoryMediator <: AbstractTheory end

"""
    make_mediator_theory(name, field_names, supertype=AbstractTheoryMediator)

Creates a type `name` representing a theory with an s-channel mediator, with
`Real`-valued fields specified by `field_names`, a comma-separated list of
symbols. The type contains a field `Γ_med` which caches the mediator width,
which is recalculated every time the other fields are updated. This field is
protected from being updated by the user with `setproperty!()`.

The supertype must be a subtype of `AbstractTheoryMediator`. Note that for some
reason `supertype()` does not work properly when called on the type created by
this macro...
"""
macro make_mediator_theory(name, field_names, supertype=AbstractTheoryMediator)
    # Get fields and constraint functions
    fields = []
    constraints = []
    for el in eval(field_names)
        if el isa Symbol
            push!(fields, :($(el)::Float64))
            push!(constraints, nothing)
        elseif length(el) == 2 && el[1] isa Symbol && eval(el[2]) isa Function
            push!(fields, :($(el[1])::Float64))
            push!(constraints, el[2])
        else
            error("'fields' entries must be :name or (:name, constraint_fn)")
        end
    end

    quote
        mutable struct $(esc(name)) <: $(esc(supertype))
            $(map(esc, fields)...)
            Γ_med::Float64  # add Γ_med field

            function $(esc(name))($(map(esc, fields)...))
                # TODO: I don't get why `[$(map(esc, fields)...)]` works, but
                # `$(map(esc, fields))` doesn't. For some reason the latter
                # gives an `Expr[]` while the former gives a `Real[]`...
                for (name, val, constr) in zip(
                    $(esc(fields)),  # needed to get field names for error
                    [$(map(esc, fields)...)],  # field values
                    [$(map(esc, constraints)...)],  # constraint functions
                )
                    if constr == nothing
                        continue
                    elseif !constr(val)
                        error("constraint '$constr' not satisfied for '$name'")
                    end
                end

                instance = new($(map(esc, fields)...), 0.0)  # set Γ_med to 0.0
                setfield!(instance, :Γ_med, Γ_med(instance))  # update Γ_med
                return instance
            end
        end

        # Prevent Γ_med from being set, and update it when fields are changed
        function Base.setproperty!(model::$(esc(name)), sym::Symbol, v)
            sym == :Γ_med && throw(error("cannot set Γ_med"))

            setfield!(model, sym, v)  # set the field if it's not width
            setfield!(model, :Γ_med, Γ_med(model))  # update width
        end
    end
end

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
    σ_tot = zero(e_cm)
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

    σ_tot = zero(e_cm)
    for fs in list_annihilation_final_states(mod)
        σ_fs = σ_χχ(e_cm, mod, fs)
        brs[fs] = σ_fs
        σ_tot += σ_fs
    end

    if σ_tot != zero(e_cm)
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
