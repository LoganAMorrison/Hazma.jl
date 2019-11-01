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
export me, ELECTRON_MASS;
export mμ, MUON_MASS;
export mπ⁰, NEUTRAL_PION_MASS;
export mπ, CHARGED_PION_MASS;
export mk⁰, NEUTRAL_KAON_MASS;
export mkL, LONG_KAON_MASS;
export mk, CHARGED_KAON_MASS;
export mη, ETA_MASS;
export mη′, ETA_PRIME_MASS;
export mρ, RHO_MASS;
export mω, OMEGA_MASS;
export mB, CHARGED_B_MASS;
export mΠ, PION_MASS_CHIRAL_LIMIT;
export mK, KAON_MASS_CHIRAL_LIMIT;
export mu, UP_QUARK_MASS;
export md, DOWN_QUARK_MASS;
export ms, STRANGE_QUARK_MASS;
export mc, CHARM_QUARK_MASS;
export mb, BOTTOM_QUARK_MASS;
export mt, TOP_QUARK_MASS;
export CM_TO_INV_MEV;
export SIGMAV_INV_MEV2_TO_CM3_PER_S;
export αem, ALPHA_EM;
export GF;
export VH;
export CMB_TEMPERATURE;
export PLANK_MASS;
export CRITICAL_ENERGY_DENSITY;
export SM_ENTROPY_DENSITY_TODAY;
export Ω_H²_CDM, OMEGA_H2_CDM;
export QU;
export QD;
export QE;
export fπ⁰, NEUTRAL_PION_DECAY_CONSTANT;
export fπ, CHARGED_PION_DECAY_CONSTANT;
export fk, CHARGED_KAON_DECAY_CONSTANT;
export B0;

include("decay/muon.jl")
export decay_spectrum_muon
include("decay/neutral_pion.jl")
export decay_spectrum_neutral_pion
include("decay/charged_pion.jl")
export decay_spectrum_charged_pion

include("positron/muon.jl")
export positron_spectrum_muon
include("positron/charged_pion.jl")
export positron_spectrum_charged_pion

include("theory.jl")
export list_annihilation_final_states, σ_χχ, br_χχ, dndeᵧ, lines_γ, dndeₑ, lines_e
include("scalar.jl")
export ScalarMediator, HiggsPortal, HeavyQuark

include("relic_density/relic_density.jl")
export solve_boltzmann
export relic_density

end # module
