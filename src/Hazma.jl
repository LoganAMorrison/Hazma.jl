module Hazma

using FastGaussQuadrature
using DelimitedFiles
using Interpolations
using SpecialFunctions
using QuadGK
using Roots
using DifferentialEquations
using ODEInterfaceDiffEq
using ForwardDiff

include("constants.jl")
export me, ELECTRON_MASS
export mμ, MUON_MASS
export mπ⁰, NEUTRAL_PION_MASS
export mπ, CHARGED_PION_MASS
export mk⁰, NEUTRAL_KAON_MASS
export mkL, LONG_KAON_MASS
export mk, CHARGED_KAON_MASS
export mη, ETA_MASS
export mη′, ETA_PRIME_MASS
export mρ, RHO_MASS
export mω, OMEGA_MASS
export mB, CHARGED_B_MASS
export mΠ, PION_MASS_CHIRAL_LIMIT
export mK, KAON_MASS_CHIRAL_LIMIT
export mu, UP_QUARK_MASS
export md, DOWN_QUARK_MASS
export ms, STRANGE_QUARK_MASS
export mc, CHARM_QUARK_MASS
export mb, BOTTOM_QUARK_MASS
export mt, TOP_QUARK_MASS
export CM_TO_INV_MEV
export SIGMAV_INV_MEV2_TO_CM3_PER_S
export αem, ALPHA_EM
export GF
export VH
export CMB_TEMPERATURE
export PLANK_MASS
export CRITICAL_ENERGY_DENSITY
export SM_ENTROPY_DENSITY_TODAY
export Ω_H²_CDM, OMEGA_H2_CDM
export QU
export QD
export QE
export fπ⁰, NEUTRAL_PION_DECAY_CONSTANT
export fπ, CHARGED_PION_DECAY_CONSTANT
export fk, CHARGED_KAON_DECAY_CONSTANT
export B0

include("utils.jl")
include("boost.jl")

include("decay/muon.jl")
export dndeᵧ_μ_decay
include("decay/neutral_pion.jl")
export dndeᵧ_π⁰_decay
include("decay/charged_pion.jl")
export decay_spectrum_charged_pion

include("positron/muon.jl")
export dndeₑ_μ_decay
include("positron/charged_pion.jl")
export dndeₑ_π_decay

include("theory.jl")
export list_annihilation_final_states, σ_χχ, br_χχ
export list_decay_final_states, Γ_med, br_med
export dndeᵧ, lines_γ
export dndeₑ, lines_e
export σ_χχ_to_ee, σ_χχ_to_μμ, σ_χχ_to_π⁰π⁰, σ_χχ_to_ππ, σ_χχ_to_γγ
# include("scalar.jl")
# export ScalarMediator, HiggsPortal, HeavyQuark

include("scalar_mediator/definitions.jl")
export ScalarMediator, HiggsPortal, HeavyQuark
include("scalar_mediator/cross_sections.jl")
export σ_χχ
export σ_χχ_to_ee
export σ_χχ_to_μμ
export σ_χχ_to_ππ
export σ_χχ_to_π⁰π⁰
export σ_χχ_to_γγ
export σ_χχ_to_ss
export σ_χχ_to_χχ
export thermal_cross_section
include("scalar_mediator/widths.jl")
export Γ_med, Γ_s_to_ee, Γ_s_to_μμ, Γ_s_to_ππ, Γ_s_to_π⁰π⁰, Γ_s_to_γγ, Γ_s_to_χχ
include("scalar_mediator/fsr.jl")
export dndeᵧ_χχ_to_eeγ, dndeᵧ_χχ_to_μμγ, dndeᵧ_χχ_to_ππγ, dndeᵧ_χχ_to_μμ
include("scalar_mediator/positron.jl")
export dndeₑ_χχ, dndeₑ_χχ_to_ππ, dndeₑ_χχ_to_μμ, dndeₑ_χχ_to_ss, lines_e

include("vector_mediator/definitions.jl")
export VectorMediator, KineticMixing
include("vector_mediator/cross_sections.jl")
export σ_χχ_to_π⁰γ, σ_χχ_to_π⁰v, σ_χχ_to_vv
include("vector_mediator/widths.jl")
export Γ_med, Γ_s_to_ee, Γ_s_to_μμ, Γ_s_to_ππ, Γ_s_to_π⁰γ, Γ_s_to_χχ


include("relic_density/relic_density.jl")
export solve_boltzmann
export relic_density

end # module
