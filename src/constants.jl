"electron mass in MeV"
const ELECTRON_MASS = 0.510998928
const me = ELECTRON_MASS
"muon mass in MeV"
const MUON_MASS = 105.6583715
const mμ = MUON_MASS
"neutral-pion mass in MeV"
const NEUTRAL_PION_MASS = 134.9766
const mπ⁰ = NEUTRAL_PION_MASS
"charged-pion mass in MeV"
const CHARGED_PION_MASS = 139.57018
const mπ = CHARGED_PION_MASS
"neutral-kaon mass in MeV"
const NEUTRAL_KAON_MASS = 497.61
const mk⁰ = NEUTRAL_KAON_MASS
"long-kaon mass in MeV"
const LONG_KAON_MASS = 497.614
const mkL = LONG_KAON_MASS
"charged kaon mass in MeV"
const CHARGED_KAON_MASS = 493.68
const mk = CHARGED_KAON_MASS
"η mass in MeV"
const ETA_MASS = 547.86
const mη = ETA_MASS
"η' mass in MeV"
const ETA_PRIME_MASS = 957.8
const mη′ = ETA_PRIME_MASS
"ρ mass in MeV"
const RHO_MASS = 775.3
const mρ = RHO_MASS
"ω mass in MeV"
const OMEGA_MASS = 782.7
const mω = OMEGA_MASS
"charged-B mass in MeV"
const CHARGED_B_MASS = 5279.29
const mB = CHARGED_B_MASS
"pion mass mass in MeV assuming chiral limit"
const PION_MASS_CHIRAL_LIMIT = (NEUTRAL_PION_MASS + CHARGED_PION_MASS) / 2
const mΠ = PION_MASS_CHIRAL_LIMIT
"kaon mass mass in MeV assuming chiral limit"
const KAON_MASS_CHIRAL_LIMIT = (NEUTRAL_KAON_MASS + CHARGED_KAON_MASS) / 2
const mK = KAON_MASS_CHIRAL_LIMIT

"up-quark mass in MeV using MS-bar"
const UP_QUARK_MASS = 2.3
const mu = UP_QUARK_MASS
"down-quark mass in MeV using MS-bar"
const DOWN_QUARK_MASS = 4.8
const md = DOWN_QUARK_MASS
"strange-quark mass in MeV using MS-bar"
const STRANGE_QUARK_MASS = 95.0
const ms = STRANGE_QUARK_MASS
"charm-quark mass in MeV using MS-bar"
const CHARM_QUARK_MASS = 1.275e3
const mc = CHARM_QUARK_MASS
"bottom-quark mass in MeV using MS-bar"
const BOTTOM_QUARK_MASS = 4.18e3
const mb = BOTTOM_QUARK_MASS
"top-quark mass in MeV using MS-bar"
const TOP_QUARK_MASS = 160.0e3
const mt = TOP_QUARK_MASS

"conversion factor from cm → 1/MeV"
const CM_TO_INV_MEV = 5.08e10
"conversion factor for converting ⟨σv⟩ in units of MeV² → cm³/s"
const SIGMAV_INV_MEV2_TO_CM3_PER_S = 1 / CM_TO_INV_MEV^2 * 3e10

"electromagnetic fine-structure constant"
const ALPHA_EM = 1.0 / 137.04
const αem = ALPHA_EM
"Fermi constant in 1/MeV²"
const GF = 1.1663787e-11
"Higgs VEV in MeV"
const VH = 246.22795e3
"CMB temperature at formation in MeV"
const CMB_TEMPERATURE = 0.235e-6
"Plank mass in MeV"
const PLANK_MASS = 1.22091e22
"critical energy density in h² MeV / cm³"
const CRITICAL_ENERGY_DENSITY = 1.05375e-2
"entropy density today 1/cm³"
const SM_ENTROPY_DENSITY_TODAY = 2891.2
"ark matter energy fraction times h²"
const OMEGA_H2_CDM = 0.1198
const Ω_H²_CDM = OMEGA_H2_CDM

"up-type-quark charge"
const QU = 2.0 / 3.0
"dow-type-quark charge"
const QD = -1.0 / 3.0
"down-type-lepton charge"
const QE = -1.0

"neutral-pion decay constant in MeV"
const NEUTRAL_PION_DECAY_CONSTANT = 91.924
const fπ⁰ = NEUTRAL_PION_DECAY_CONSTANT
"charged-pion decay constant in MeV"
const CHARGED_PION_DECAY_CONSTANT = 92.2138
const fπ = CHARGED_PION_DECAY_CONSTANT
"charged-kaon decay constant in MeV"
const CHARGED_KAON_DECAY_CONSTANT = 110.379
const fk = CHARGED_KAON_DECAY_CONSTANT
"quark-condinsate constant: B₀³ = -⟨q̄q⟩/(3F₀²)"
const B0 = PION_MASS_CHIRAL_LIMIT^2 / (UP_QUARK_MASS + DOWN_QUARK_MASS)
"PDG convension of charged-pion decay constant in MeV"
const CHARGED_PION_DECAY_CONSTANT_PDG = 130.2
const fπ_PDG = CHARGED_PION_DECAY_CONSTANT_PDG
"PDG convension of charged-kaon decay constant in MeV"
const CHARGED_KAON_DECAY_CONSTANT_PDG = 155.6
const fk_PDG = CHARGED_KAON_DECAY_CONSTANT_PDG

# const G8 = 5.47
# const G27 = 0.392
# const gv = 67.0
# const fv = 153.0

"""
INTERNAL USE ONLY
"""

"up-down CMK"
const VUD = 0.974267
"up-strange CMK"
const VUS = 0.2248
"top-strange CMK"
const VTS = -0.0405 - 0.00075987im
"top-bottom CMK"
const VTB = 0.999139
"top-down CMK"
const VTD = 0.00823123 - 0.00328487im

"decay-width of ρ(770) in MeV"
const WIDTH_RHO = 149.1
const Γρ = WIDTH_RHO
"decay-width of charged-B in MeV"
const WIDTH_B = 4.018e-10
const ΓB = WIDTH_B
"decay-width of neutral-B in MeV"
const WIDTH_B0 = 4.333e-10
const ΓB0 = WIDTH_B0
"decay-width of charged kaon in MeV"
const WIDTH_K = 5.32e-14
const Γk = WIDTH_K
"decay-width of long-kaon in MeV"
const WIDTH_KL = 1.29e-14
const ΓkL = WIDTH_KL
"decay-width of short-kaon in MeV "
const WIDTH_KS = 7.3510e-12
const ΓkS = WIDTH_KS
"decay-width of charged-pion in MeV"
const WIDTH_PI = 2.5284e-14
const Γπ = WIDTH_PI
"decay-width of neutral-pion in MeV"
const WIDTH_PI0 = 7.73e-6
const Γπ⁰ = WIDTH_PI0
"decay-width of η in MeV"
const WIDTH_ETA = 1.31e-3
const Γη = WIDTH_ETA
"decay-width of η'(958) in MeV"
const WIDTH_ETAPRIME = 0.196
const Γη = WIDTH_ETA

"axial-vector form-factor of the charged-pion"
const AXIAL_FORM_FACTOR_PI = 0.0119
"vector form-factor of the charged-pion"
const VECTOR_FORM_FACTOR_PI = 0.0254
"slope of energy dependence of Vπ(s) = Vπ(0)(1 + a⋅s), s = 1 - 2Eγ/mπ"
const VECTOR_FORM_FACTOR_SLOPE_PI = 0.095
"axial-vector form-factor of the charged-kaon"
const AXIAL_FORM_FACTOR_KAON = 0.042
"vector form-factor of the charged-kaon"
const VECTOR_FORM_FACTOR_KAON = 0.096

"""Branching Fractions"""

"branching-ratio π⁰ → γγ"
const BR_PI0_TO_GG = 0.9882
const br_π⁰_γγ = BR_PI0_TO_GG

"branching-ratio π⁻ → μ⁻ν"
const BR_PI_TO_MUNU = 0.9998
const br_π_μν = BR_PI_TO_MUNU
"branching-ratio π⁻ → e⁻ν"
const BR_PI_TO_ENU = 0.000123
const br_π_eν = BR_PI_TO_ENU

"branching-ratio Ks → ππ"
const BR_KS_TO_PIPI = 0.6920
const br_kS_ππ = BR_KS_TO_PIPI
"branching-ratio Ks → π⁰π⁰"
const BR_KS_TO_PI0PI0 = 0.3069
const br_kS_π⁰π⁰ = BR_KS_TO_PI0PI0

"branching-ratio KL → πeν"
const BR_KL_TO_PIENU = 0.4055
const br_kL_πeν = BR_KL_TO_PIENU
"branching-ratio KL → πμν"
const BR_KL_TO_PIMUNU = 0.2704
const br_kL_πμν = BR_KL_TO_PIMUNU
"branching-ratio KL → π⁰π⁰π⁰"
const BR_KL_TO_3PI0 = 0.1952
const br_kL_3π⁰ = BR_KL_TO_3PI0
"branching-ratio KL → πππ⁰"
const BR_KL_TO_2PIPI0 = 0.1254
const br_kLπππ⁰ = BR_KL_TO_2PIPI0

"branching-ratio K⁻ → μ⁻ν"
const BR_K_TO_MUNU = 0.6356
const br_k_μν = BR_K_TO_MUNU
"branching-ratio K⁻ → π⁻π⁰"
const BR_K_TO_PIPI0 = 0.2067
const br_k_ππ⁰ = BR_K_TO_PIPI0
"branching-ratio K⁻ → π⁻π⁺π⁻"
const BR_K_TO_3PI = 0.05583
const br_k_3π = BR_K_TO_3PI
"branching-ratio K⁻ → π⁰e⁻ν"
const BR_K_TO_PI0ENU = 0.0507
const br_k_π⁰eν = BR_K_TO_PI0ENU
"branching-ratio K⁻ → π⁰μ⁻ν"
const BR_K_TO_PI0MUNU = 0.03352
const br_k_π⁰μν = BR_K_TO_PI0MUNU
"branching-ratio K⁻ → π⁻π⁰π⁰"
const BR_K_TO_2PIPI0 = 0.01760
const br_kπππ⁰ = BR_K_TO_2PIPI0

"branching-ratio η → γγ"
const BR_ETA_TO_GG = 0.3941
const br_η_γγ = BR_ETA_TO_GG
"branching-ratio η → π⁰π⁰π⁰"
const BR_ETA_TO_3PI0 = 0.3268
const br_η_3π⁰ = BR_ETA_TO_3PI0
"branching-ratio η → π⁻π⁺π⁰"
const BR_ETA_TO_2PIPI0 = 0.2292
const br_ηπππ⁰ = BR_ETA_TO_2PIPI0
"branching-ratio η → π⁺π⁻γ"
const BR_ETA_TO_2PIG = 0.0422
const br_ηππγ = BR_ETA_TO_2PIG

"branching-ratio η' → π⁺π⁻η"
const BR_ETAP_TO_2PIETA = 0.429
const br_η′ππη = BR_ETAP_TO_2PIETA
"branching-ratio η' → ργ"
const BR_ETAP_TO_RHOG = 0.291
const br_η′ππη = BR_ETAP_TO_2PIETA
"branching-ratio η' → π⁰π⁰η"
const BR_ETAP_2PI0ETA = 0.222
const br_η′_π⁰π⁰η = BR_ETAP_2PI0ETA
"branching-ratio η' → ωγ"
const BR_ETAP_TO_OMEGAG = 0.0275
const br_η′_ωγ = BR_ETAP_TO_OMEGAG
"branching-ratio η' → γγ"
const BR_ETAP_TO_GG = 0.0220
const br_η′_γγ = BR_ETAP_TO_GG
"branching-ratio η' → π⁰π⁰π⁰"
const BR_ETAP_TO_3PI0 = 0.0214
const br_η′_3π⁰ = BR_ETAP_TO_3PI0
"branching-ratio η' → μ⁻μ⁺γ"
const BR_ETAP_TO_MUMUG = 0.0108
const br_η′_μμγ = BR_ETAP_TO_MUMUG

"branching-ratio ω → π⁻π⁺π⁰"
const BR_OMEGA_TO_2PIPI0 = 0.892
const br_ωπππ⁰ = BR_OMEGA_TO_2PIPI0
"branching-ratio ω → π⁰γ"
const BR_OMEGA_TO_PI0G = 0.0828
const br_ω_π⁰γ = BR_OMEGA_TO_PI0G
"branching-ratio ω → π⁻π⁺"
const BR_OMEGA_TO_2PI = 0.0153
const br_ωππ = BR_OMEGA_TO_2PI
