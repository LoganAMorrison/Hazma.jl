"electron mass in MeV"
const ELECTRON_MASS = 0.510998928;
"muon mass in MeV"
const MUON_MASS = 105.6583715;
"neutral-pion mass in MeV"
const NEUTRAL_PION_MASS = 134.9766;
"charged-pion mass in MeV"
const CHARGED_PION_MASS = 139.57018;
"neutral-kaon mass in MeV"
const NEUTRAL_KAON_MASS = 497.61;
"long-kaon mass in MeV"
const LONG_KAON_MASS = 497.614;
"charged kaon mass in MeV"
const CHARGED_KAON_MASS = 493.68;
"η mass in MeV"
const ETA_MASS = 547.86;
"η' mass in MeV"
const ETA_PRIME_MASS = 957.8;
"ρ mass in MeV"
const RHO_MASS = 775.3;
"ω mass in MeV"
const OMEGA_MASS = 782.7;
"charged-B mass in MeV"
const CHARGED_B_MASS = 5279.29;
"pion mass mass in MeV assuming chiral limit"
const PION_MASS_CHIRAL_LIMIT = (NEUTRAL_PION_MASS + CHARGED_PION_MASS) / 2;
"kaon mass mass in MeV assuming chiral limit"
const KAON_MASS_CHIRAL_LIMIT = (NEUTRAL_KAON_MASS + CHARGED_KAON_MASS) / 2;

"up-quark mass in MeV using MS-bar"
const UP_QUARK_MASS = 2.3;
"down-quark mass in MeV using MS-bar"
const DOWN_QUARK_MASS = 4.8;
"strange-quark mass in MeV using MS-bar"
const STRANGE_QUARK_MASS = 95.0;
"charm-quark mass in MeV using MS-bar"
const CHARM_QUARK_MASS = 1.275e3;
"bottom-quark mass in MeV using MS-bar"
const BOTTOM_QUARK_MASS = 4.18e3;
"top-quark mass in MeV using MS-bar"
const TOP_QUARK_MASS = 160.0e3;

"conversion factor from cm → 1/MeV"
const CM_TO_INV_MEV = 5.08e10;
"conversion factor for converting ⟨σv⟩ in units of MeV² → cm³/s"
const SIGMAV_INV_MEV2_TO_CM3_PER_S = 1 / CM_TO_INV_MEV^2 * 3e10;

"electromagnetic fine-structure constant"
const ALPHA_EM = 1.0 / 137.04;
"Fermi constant in 1/MeV²"
const GF = 1.1663787e-11;
"Higgs VEV in MeV"
const VH = 246.22795e3;
"CMB temperature at formation in MeV"
const CMB_TEMPERATURE = 0.235e-6;
"Plank mass in MeV"
const PLANK_MASS = 1.22091e22;
"critical energy density in h² MeV / cm³"
const CRITICAL_ENERGY_DENSITY = 1.05375e-2;
"entropy density today 1/cm³"
const SM_ENTROPY_DENSITY_TODAY = 2891.2;
"ark matter energy fraction times h²"
const OMEGA_H2_CDM = 0.1198;

"up-type-quark charge"
const QU = 2.0 / 3.0;
"dow-type-quark charge"
const QD = -1.0 / 3.0;
"down-type-lepton charge"
const QE = -1.0;

"neutral-pion decay constant in MeV"
const NEUTRAL_PION_DECAY_CONSTANT = 91.924;
"charged-pion decay constant in MeV"
const CHARGED_PION_DECAY_CONSTANT = 92.2138;
"charged-kaon decay constant in MeV"
const CHARGED_KAON_DECAY_CONSTANT = 110.379;
"quark-condinsate constant: B₀³ = -⟨q̄q⟩/(3F₀²)"
const B0 = PION_MASS_CHIRAL_LIMIT^2 / (UP_QUARK_MASS + DOWN_QUARK_MASS);
"PDG convension of charged-pion decay constant in MeV"
const CHARGED_PION_DECAY_CONSTANT_PDG = 130.2;
"PDG convension of charged-kaon decay constant in MeV"
const CHARGED_KAON_DECAY_CONSTANT_PDG = 155.6;

# const G8 = 5.47
# const G27 = 0.392
# const gv = 67.0
# const fv = 153.0

"""
INTERNAL USE ONLY
"""

"up-down CMK"
const VUD = 0.974267;
"up-strange CMK"
const VUS = 0.2248;
"top-strange CMK"
const VTS = -0.0405 - 0.00075987im;
"top-bottom CMK"
const VTB = 0.999139;
"top-down CMK"
const VTD = 0.00823123 - 0.00328487im;

"decay-width of ρ(770) in MeV"
const WIDTH_RHO = 149.1;
"decay-width of charged-B in MeV"
const WIDTH_B = 4.018e-10;
"decay-width of neutral-B in MeV"
const WIDTH_B0 = 4.333e-10;
"decay-width of charged kaon in MeV"
const WIDTH_K = 5.32e-14;
"decay-width of long-kaon in MeV"
const WIDTH_KL = 1.29e-14;
"decay-width of short-kaon in MeV "
const WIDTH_KS = 7.3510e-12;
"decay-width of charged-pion in MeV"
const WIDTH_PI = 2.5284e-14;
"decay-width of neutral-pion in MeV"
const WIDTH_PI0 = 7.73e-6;
"decay-width of η in MeV"
const WIDTH_ETA = 1.31e-3;
"decay-width of η'(958) in MeV"
const WIDTH_ETAPRIME = 0.196;

"axial-vector form-factor of the charged-pion"
const AXIAL_FORM_FACTOR_PI = 0.0119;
"vector form-factor of the charged-pion"
const VECTOR_FORM_FACTOR_PI = 0.0254;
"slope of energy dependence of Vπ(s) = Vπ(0)(1 + a⋅s), s = 1 - 2Eγ/mπ"
const VECTOR_FORM_FACTOR_SLOPE_PI = 0.095;
"axial-vector form-factor of the charged-kaon"
const AXIAL_FORM_FACTOR_KAON = 0.042;
"vector form-factor of the charged-kaon"
const VECTOR_FORM_FACTOR_KAON = 0.096;


"""Branching Fractions"""


"branching-ratio π⁰ → γγ"
const BR_PI0_TO_GG = 0.9882;

"branching-ratio π⁻ → μ⁻ν"
const BR_PI_TO_MUNU = 0.9998;
"branching-ratio π⁻ → e⁻ν"
const BR_PI_TO_ENU = 0.000123;

"branching-ratio Ks → ππ"
const BR_KS_TO_PIPI = 0.6920;
"branching-ratio Ks → π⁰π⁰"
const BR_KS_TO_PI0PI0 = 0.3069;

"branching-ratio KL → πeν"
const BR_KL_TO_PIENU = 0.4055;
"branching-ratio KL → πμν"
const BR_KL_TO_PIMUNU = 0.2704;
"branching-ratio KL → π⁰π⁰π⁰"
const BR_KL_TO_3PI0 = 0.1952;
"branching-ratio KL → πππ⁰"
const BR_KL_TO_2PIPI0 = 0.1254;

"branching-ratio K⁻ → μ⁻ν"
const BR_K_TO_MUNU = 0.6356;
"branching-ratio K⁻ → π⁻π⁰"
const BR_K_TO_PIPI0 = 0.2067;
"branching-ratio K⁻ → π⁻π⁺π⁻"
const BR_K_TO_3PI = 0.05583;
"branching-ratio K⁻ → π⁰e⁻ν"
const BR_K_TO_PI0ENU = 0.0507;
"branching-ratio K⁻ → π⁰μ⁻ν"
const BR_K_TO_PI0MUNU = 0.03352;
"branching-ratio K⁻ → π⁻π⁰π⁰"
const BR_K_TO_PI2PI0 = 0.01760;

"branching-ratio η → γγ"
const BR_ETA_TO_GG = 0.3941;
"branching-ratio η → π⁰π⁰π⁰"
const BR_ETA_TO_3PI0 = 0.3268;
"branching-ratio η → π⁻π⁺π⁰"
const BR_ETA_TO_2PIPI0 = 0.2292;
"branching-ratio η → π⁺π⁻γ"
const BR_ETA_TO_2PIG = 0.0422;

"branching-ratio η' → π⁺π⁻η"
const BR_ETAP_TO_2PIETA = 0.429;
"branching-ratio η' → ργ"
const BR_ETAP_TO_RHOG = 0.291;
"branching-ratio η' → π⁰π⁰η"
const BR_ETAP_2PI0ETA = 0.222;
"branching-ratio η' → ωγ"
const BR_ETAP_TO_OMEGAG = 0.0275;
"branching-ratio η' → γγ"
const BR_ETAP_TO_GG = 0.0220;
"branching-ratio η' → π⁰π⁰π⁰"
const BR_ETAP_TO_3PI0 = 0.0214;
"branching-ratio η' → μ⁻μ⁺γ"
const BR_ETAP_TO_MUMUG = 0.0108;

"branching-ratio ω → π⁻π⁺π⁰"
const BR_OMEGA_TO_2PIPI0 = 0.892;
"branching-ratio ω → π⁰γ"
const BR_OMEGA_TO_PI0G = 0.0828;
"branching-ratio ω → π⁻π⁺"
const BR_OMEGA_TO_2PI = 0.0153;
