from types import SimpleNamespace

con = SimpleNamespace()

# Ion concentrations
con.Mgi = 2.5          # [mM] intracellular Mg+
con.Nao = 140          # [mM] extracellular Na+
con.Cao = 2            # [mM] extracellular Ca2+
con.Ko = 5.4           # [mM] extracellular K+
con.Nai = 10           # [mM] intracellular Na+
con.Ki = 139.8854603066  # [mM] intracellular K+

# Cell compartments
con.C = 0.025          # [nF] cell electric capacitance
con.V_cell = 3         # [pL]
con.V_sub = 0.03328117 # [pL] submembrane cell volume (0.035097874)
con.V_jSR = 0.0036
con.V_i = 1.34671883   # [pL] myoplasmic volume
con.V_nSR = 0.0348     # [pL] volume of network SR (Ca2+ uptake store)
con.F = 96485          # C/M is the Faraday constant
con.T = 310.5          # K is the absolute temperature for 37°C
con.R = 8.314472       # J/(M K) is the universal gas constant
con.E_st = 17          # [mV] apparent reversal potential of I_st
con.E_CaL = 47         # [mV] apparent reversal potential of I_CaL
con.ENa15 = 41.5761    # [mV] Reversal potential of I_Na1.5
con.E_CaT = 45         # [mV] apparent reversal potential of I_CaT
con.E_T = 1000 * (con.R * con.T / con.F)  # [mV] 'RT/F' factor, 1000 converts V to mV

# Na+/Ca2+ exchanger current parameters
con.K_1ni = 395.3      # [mM] Nai binding to the first site on NCX
con.K_1no = 1628       # [mM] Nao binding to first site on NCX
con.K_2ni = 2.289      # [mM] Nai binding to second site on NCX
con.K_2no = 561.4      # [mM] Nao binding to second site on NCX
con.K_3ni = 26.44      # [mM] Nai binding to third site on NCX
con.K_3no = 4.663      # [mM] Nao binding to third site on NCX
con.K_ci = 0.0207      # [mM] Cai binding on NCX transporter
con.K_co = 3.663       # [mM] Cao binding to NCX transporter
con.K_cni = 26.44      # [mM] Nai and Cai simultaneous binding to NCX
con.Q_ci = 0.1369      # Cai occlusion reaction of NCX
con.Q_co = 0           # Cao occlusion reaction of NCX
con.Q_n = 0.4315       # Na occlusion of NCX

# Ca2+ flux parameters
con.tho_difCa = 0.04   # [ms] Time constant of Ca diffusion from the submembrane to myoplasm
con.TC_tot = 0.031     # [mM] Total concentration of the troponin-Ca site
con.TMC_tot = 0.062    # [mM] Total concentration of the troponin-Mg site
con.k_fTC = 88.8       # [mM/ms] Ca association constant for troponin
con.k_fTMC = 237.7     # [mM/ms] Ca association constant for the troponin-Mg site
con.k_bTC = 0.446      # [1/ms] Ca dissociation constant for the troponin-Ca site
con.k_bTMC = 0.00751   # [1/ms] Ca dissociation constant for the troponin-Mg site
con.k_fTMM = 2.277     # [mM/ms] Mg association constant for the troponin-Mg site
con.k_bTMM = 0.751     # [1/ms] Mg dissociation constant for the troponin-Mg site
con.CM_tot = 0.045     # [mM] Total calmodulin concentration
con.k_fCM = 227.7      # [1/(mM*ms)] Ca association constant for calmodulin 
con.k_bCM = 0.542      # [1/ms] Ca dissociation constant for calmodulin
con.k_bCQ = 0.445      # [1/ms] Ca dissociation constant for calsequestrin
con.k_iCa = 0.5        # [1/(mM*ms)]
con.k_im = 0.005       # [1/ms]
con.k_oCa_max = 10     # [1/(mM^2*ms)]
con.RyR_min = 0.0127   # derived from PP1 activity
con.RyR_max = 0.02
con.n_RyR = 9.773
con.k_05Ry = 0.7
con.EC_50SR = 0.45     # [mM]
con.MaxSR = 15
con.MinSR = 1
con.HSR = 2.5
con.n_up = 2
con.P_up = 0.02        # [mM/ms] Rate constant for Ca uptake by the Ca pump in the network SR
con.k_s = 250e3        # [ms^-1]
con.K_mf = 0.00008     # [mM]
con.K_mr = 4.5         # [mM]
con.tho_tr = 40        # [ms] Time constant for Ca transport from the network to junctional SR
con.K_up = 0.6e-3      # [mM] Half-maximal Cai for Ca uptake in the network SR
con.K_mfCa = 0.00035   # [mM] Dissociation constant of Ca2+-dependent I_CaL inactivation
con.alpha_fCa = 0.021  # [1/ms] Ca2+ dissociation rate constant for I_CaL

# Membrane parameters
con.g_sus = 0.00039060               # [nS/pF]*[nF]=[µS] Normalized conductance for I_sus channels
con.g_st = 0.000000006               # [nS/pF]*[nF]=[µS] Normalized conductance for I_st channels
con.g_bNa = 0.0001215                # [nS/pF]*[nF]=[µS]
con.g_K1 = 0.229 * 0.0039228 * 0.9   # [nS/pF]*[nF]=[µS]
con.g_Ks = 0.000299                  # [nS/pF]*[nF]=[µS] Normalized conductance for I_Ks channels
con.g_Na15 = 0.000237 * 0.025        # [nS/pF]*[nF]=[µS]
con.g_Na11 = 0.000237 * 0.025        # [nS/pF]*[nF]=[µS]
con.g_CaL12 = 0.2832 * 0.025         # [nS/pF]*[nF]=[µS]
con.g_CaL13 = 0.9936 * 0.025         # [nS/pF]*[nF]=[µS]
con.g_CaT = 0.5600 * 0.025           # [nS/pF]*[nF]=[µS] Normalized conductance for I_CaT channels
con.g_If = 0.228 * 0.025             # [nS/pF]*[nF]=[µS]
con.g_Kr = 0.0960 * 0.025            # [nS/pF]*[nF]=[µS]
con.g_to = 0.00492                   # [nS/pF]*[nF]=[µS] Normalized conductance for I_to channels
con.K_NaCa = 220 * 0.025             # [pA/pF]
con.g_bCa = 0.0006 * 0.025           # [nS/pF]*[nF]=[µS]
con.g_bK = 0.0001 * 0.025            # [nS/pF]*[nF]=[µS]
con.I_NaKmax = 1.85 * 0.077          # [pA/pF] Maximal Na/K pump current conductance
con.K_mNa = 14                       # [mM] half-maximal Nai for I_NaK
con.K_mK = 1.4                       # [mM] half-maximal Ko for I_NaK
con.g_KACh_max = 2 * 0.14241818

# AC–cAMP/PKA signaling parameters within the numerical model
con.K_ACI = 0.016        # [min^-1] Non-Ca2+ AC activity
con.K_AC = 0.0735        # [min^-1] Non-Ca2+ AC activity
con.K_Ca = 0.000178      # [mM] Maximal Ca2+ AC activation
con.K_AC_Ca = 0.000024   # [mM] Half-maximal Ca2+ AC activation
con.k_PKA = 9000         # [pmol/protein/min] Maximal PKA activity
con.k_PKA_cAMP = 284.5   # [pmol/protein] Half-maximal PKA activation
con.n_PKA = 5            # Hill coefficient
con.PKA_PLB = 1          # New - Limor (for PKA-PLB inhibition)

# Phosphorylation parameters
con.k_PLBp = 52.25       # [1/min] Maximal PLB phosphorylation
con.n_PLB = 1            # Hill coefficient
con.k_PKA_PLB = 1.651    # Half-maximal PLB phosphorylation
con.PP1 = 0.89           # [µM] PP1 concentration
con.k_PP1 = 23.575       # [1/µM/min] Maximal PP1 activity
con.k_pp1_PLB = 0.06967  # Half maximal PP1 activity

# Force parameters
con.SL_lo = 0.8e-6       # [m] A constant coefficient that describes the effect of the actin- and myosin-filament lengths on the single overlap length. 
con.Nc = 2e13            # [1/mm^2] SAN cross-section area
con.FK0 = 350            # [1/mM] Cross-bridge independent coefficient of calcium affinity
con.FN = 3.5             # Hill coefficient
con.FK05 = 2.5e9         # [1/mm^3] Half-maximal cross-bridge Ca2+ affinity
con.Fkl = 60             # [1/mM/ms] Rate constant of Ca binding to troponin low-affinity sites
con.Ff = 40e-3           # [1/ms] The cross-bridge turnover rate from weak to strong conformation
con.Fg0 = 30e-3          # [1/ms] The cross-bridge weakening rate at isometric regime
con.Fg1 = 4.4e6          # [1/m] Mechanical-feedback coefficient
con.Fxb = 2e-9           # [mN] Unitary force per cross-bridge at isometric regime

# ATP parameters
con.ATP_max = 2.533      # [mM]
con.kATP = 61.42 * 100   # []
con.k_ATP05 = 6724       # []
con.cAMPb = 20           # [pmol/mg protein] Baseline cAMP
con.K_ATPmin = 6034      # []
con.n_ATP = 3.36         # []

# Perturbations
con.ISO = 0
con.CCh = 0
con.IBMX = 0
con.ISO_power = 1.5
con.PKA_temp = 3.97
