module Hazma

using FastGaussQuadrature

include("constants.jl")
export ELECTRON_MASS;
export MUON_MASS5;
export NEUTRAL_PION_MASS;
export CHARGED_PION_MASS;
export NEUTRAL_KAON_MASS;
export LONG_KAON_MASS;
export CHARGED_KAON_MASS;
export ETA_MASS;
export ETA_PRIME_MASS;
export RHO_MASS;
export OMEGA_MASS;
export CHARGED_B_MASS;
export PION_MASS_CHIRAL_LIMIT;
export KAON_MASS_CHIRAL_LIMIT;
export UP_QUARK_MASS;
export DOWN_QUARK_MASS;
export STRANGE_QUARK_MASS;
export CHARM_QUARK_MASS;
export BOTTOM_QUARK_MASS;
export TOP_QUARK_MASS;
export CM_TO_INV_MEV;
export SIGMAV_INV_MEV2_TO_CM3_PER_S;
export ALPHA_EM;
export GF;
export VH;
export CMB_TEMPERATURE;
export PLANK_MASS;
export CRITICAL_ENERGY_DENSITY;
export SM_ENTROPY_DENSITY_TODAY;
export OMEGA_H2_CDM;
export QU;
export QD;
export QE;
export NEUTRAL_PION_DECAY_CONSTANT;
export CHARGED_PION_DECAY_CONSTANT;
export CHARGED_KAON_DECAY_CONSTANT;
export B0;

include("decay/muon.jl")
export decay_spectrum_muon
include("decay/neutral_pion.jl")
export decay_spectrum_neutral_pion
include("decay/charged_pion.jl")
export decay_spectrum_charged_pion


end # module
