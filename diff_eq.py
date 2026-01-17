import numpy as np
from constants import con         
from update_PKA import update_pka

def diff_eq(t, Y):

    Y = np.asarray(Y, dtype=float)

    Vm  = Y[0]   # [mV] Membrane potential
    qa  = Y[1]   # I_st activation
    qi  = Y[2]   # I_st inactivation
    dT  = Y[3]   # I_CaT inactivation gate
    fT  = Y[4]   # I_CaT activation gate
    pa  = Y[5]   # I_Kr activation
    pi  = Y[6]   # I_Kr inactivation
    xs  = Y[7]   # I_Ks activation
    fL12 = Y[8]  # I_CaL1.2 inactivation
    dL12 = Y[9]  # I_CaL1.2 activation
    fL13 = Y[10] # I_CaL1.3 inactivation
    dL13 = Y[11] # I_CaL1.3 activation
    fCa  = Y[12] # Ca2+-dependent inactivation gating variable for I_(CaL,1.2)and I_(CaL,1.3)
    r    = Y[13] # Activation gating variable of I_toand I_sus
    q    = Y[14] # Inactivation gating variable of I_to
    m15  = Y[15] # Activation gating variable of Nav1.5
    h15  = Y[16] # Fast inactivation gating variable of Nav1.5
    j15  = Y[17] # Slow inactivation gating variable of Nav1.5
    m11  = Y[18] # Activation gating variable of Nav1.1
    h11  = Y[19] # Fast inactivation gating variable of Nav1.1
    j11  = Y[20] # Slow inactivation gating variable of Nav1.1
    y    = Y[21] # Activation gating variable of I_f
    Cai  = Y[22] # [mM] Intracellular Ca2+ concentration or Ca2+ concentration in the cytosol
    CajSR = Y[23] # [mM] Ca2+ concentration in the JSR
    CanSR = Y[24] # [mM] Ca2+ concentration in the NSR
    Casub = Y[25] # [mM] Ca2+ concentration in the subspace
    f_TC  = Y[26] # Fractional occupancy of the troponin Ca2+ site by [Ca2+]i
    f_TMC = Y[27] # Fractional occupancy of the troponin Mg2+ site by [Ca2+]i
    f_TMM = Y[28] # Fractional occupancy of the troponin Mg2+ site by Mg2+
    f_CMs = Y[29] # Fractional occupancy of calmodulin by [Ca2+]sub
    f_CMi = Y[30] # Fractional occupancy of calmodulin by [Ca2+]i
    f_CQ  = Y[31] # Fractional occupancy of calsequestrin by [Ca2+]rel
    R     = Y[32] # RyR reactivated
    OO    = Y[33] # RyR open
    S     = Y[34] # RyR inactive
    RI    = Y[35] # RyR inactivated
    w     = Y[36] # IKACh gate
    cAMP  = Y[37] # cAMP
    PLB   = Y[38] # PLB phosphorylation level
    A     = Y[39] # The density of regulatory units with bound Ca2+ and adjacent weak cross-bridges
    TT    = Y[40] # The density of regulatory units with bound Ca2+ and adjacent strong cross-bridge
    U     = Y[41] # The density of regulatory units without bound but with adjacent strong cross-bridge
    SL    = Y[42] # Sarcomere length


    # ============== S.3.1 General ==============
    # S.3.1 Reversal potentials 
    ENa = con.E_T * np.log(con.Nao / con.Nai)
    EK  = con.E_T * np.log(con.Ko / con.Ki)
    EKs = con.E_T * np.log((con.Ko + 0.12 * con.Nao) / (con.Ki + 0.12 * con.Nai))
    ECa = (con.E_T / 2.0) * np.log(con.Cao / Casub)

    # ============== S.3.2 AC–cAMP–PKA signaling ==============
    # S.3.2.1 PKA activity
    PKA = update_pka(cAMP)

    # S.3.2.2 ATP–ADP
    ATP = (con.ATP_max *((con.kATP * (cAMP * 100.0 / con.cAMPb) ** con.n_ATP)/ (con.k_ATP05 + (cAMP * 100.0 / con.cAMPb) ** con.n_ATP)- con.K_ATPmin)/ 100.0)

    # S.3.2.3 cAMP dynamics
    k_iso = 0.1599 * (con.ISO ** 1.5 / (76.5441 ** 0.6238 + con.ISO ** 1.5))
    k_ibmx = 1.0-(0.86*con.IBMX)/(con.IBMX-3.46)
    k_CCh = 0.0146 * (con.CCh ** 1.4402 / (51.7331 ** 1.4402 + con.CCh ** 1.4402))
    k_1 = con.K_ACI + con.K_AC / (1.0 + np.exp((con.K_Ca - con.k_bCM * f_CMi / (con.k_fCM * (1.0 - f_CMi))) / con.K_AC_Ca))
    k_2 = k_ibmx * 265.3512 * (cAMP ** 5.7343) / (24.7290 ** 6.7343 + cAMP ** 6.7343)
    k_3 = (con.k_PKA * cAMP ** (con.n_PKA - 1.0)) / (con.k_PKA_cAMP ** con.n_PKA + cAMP ** con.n_PKA)
    dcAMP = (k_iso * (ATP * 0.6 * 1e3)+ k_1 * (ATP * 0.6 * 1e3)- k_2 * cAMP- k_3 * cAMP - k_CCh * (ATP * 0.6 * 1e3)) / 60000.0

    # S.3.2.4 PLB activity
    k_4 = ((con.k_PLBp * (con.PKA_PLB * PKA) ** con.n_PLB)/ (con.k_PKA_PLB ** con.n_PLB + (con.PKA_PLB * PKA) ** con.n_PLB))
    k_5 = con.k_PP1 * con.PP1 * PLB / (con.k_pp1_PLB + PLB)
    dPLB = (k_4 - k_5) / 60000.0

    # ============== S.3.3 Membrane currents ==============

    # S.3.3.1 4-aminopyridine-sensitive currents, I_toand I_sus
    It0  = con.g_to  * (Vm - EK) * q * r
    Isus = con.g_sus * (Vm - EK) * r
    q_inf = 1.0 / (1.0 + np.exp((Vm + 49.0) / 13.0))
    r_inf = 1.0 / (1.0 + np.exp(-(Vm - 19.3) / 15.0))
    tau_q = (6.06 + 39.102 / (0.57 * np.exp(-0.08 * (Vm + 44.0)) + 0.065 * np.exp(0.1 * (Vm + 45.93)))) / 0.67
    tau_r = (2.75 + 14.40516 / (1.037 * np.exp(0.09 * (Vm + 30.61)) + 0.369 * np.exp(-0.12 * (Vm + 23.84)))) / 0.303
    dr = (r_inf - r) / tau_r
    dq = (q_inf - q) / tau_q

    # S.3.3.2 Ca2+ background current, IbCa
    Ib_Ca = con.g_bCa * (Vm - ECa)
    
    # S.3.3.3 K+ background current, IbK
    Ib_K  = con.g_bK  * (Vm - EK)
    
    # S.3.3.4 Na+ background current, IbNa
    Ib_Na = con.g_bNa * (Vm - ENa)

    # S.3.3.5 L-type channel current, ICaL
    b_CaL = -0.2152 + 1.6913 * PKA ** 10.0808 / (0.8836 ** 10.0808 + PKA ** 10.0808)
    ICaL12 = con.g_CaL12 * (1.0 + b_CaL) * (Vm - con.E_CaL) * dL12 * fL12 * fCa
    ICaL13 = con.g_CaL13 * (1.0 + b_CaL) * (Vm - con.E_CaL) * dL13 * fL13 * fCa
    ICaL   = ICaL12 + ICaL13

    dL12_inf = 1.0 / (1.0 + np.exp(-(Vm + 3.0)  / 5.0))
    fL12_inf = 1.0 / (1.0 + np.exp((Vm + 36.0) / 4.6))
    dL13_inf = 1.0 / (1.0 + np.exp(-(Vm + 13.5) / 6.0))
    fL13_inf = 1.0 / (1.0 + np.exp((Vm + 35.0) / 7.3))
    fCa_inf = con.K_mfCa / (con.K_mfCa + Casub)

    alpha_dL = -28.39 * (Vm + 35.0) / (np.exp(-(Vm + 35.0) / 2.5) - 1.0) + 408.173
    beta_dL  = 11.43  * (Vm - 5.0)  / (np.exp(0.4 * (Vm - 5.0)) - 1.0)
    
    tau_dL = 2000.0 / (alpha_dL + beta_dL)
    tau_fL = 7.4 + 45.77 * np.exp(-0.5 * (Vm + 28.1) ** 2 / (11.0 * 11.0))
    tau_fCa = fCa_inf / con.alpha_fCa

    dfL12 = (fL12_inf - fL12) / tau_fL
    ddL12 = (dL12_inf - dL12) / tau_dL
    dfL13 = (fL13_inf - fL13) / tau_fL
    ddL13 = (dL13_inf - dL13) / tau_dL
    dfCa  = (fCa_inf  - fCa)  / tau_fCa

    # S.3.3.6 T-type Ca2+ current, ICaT
    ICaT = con.g_CaT * (Vm - con.E_CaT) * dT * fT
    dT_inf = 1.0 / (1.0 + np.exp(-(Vm + 26.0) / 6.0))
    fT_inf = 1.0 / (1.0 + np.exp((Vm + 61.7) / 5.6))
    tau_dT = 1.0 / (1.068 * np.exp((Vm + 26.3) / 30.0) + 1.068 * np.exp(-(Vm + 26.3) / 30.0))
    tau_fT = 1.0 / (0.0153 * np.exp(-(Vm + 61.7) / 83.3) + 0.015  * np.exp((Vm + 61.7) / 15.38))
    ddT = (dT_inf - dT) / tau_dT
    dfT = (fT_inf - fT) / tau_fT

    # S.3.3.7 Hyperpolarization-activated, funny current, If
    IfNa = 0.3833 * con.g_If * (Vm - ENa) * y
    IfK  = 0.6167 * con.g_If * (Vm - EK)  * y
    If   = IfNa + IfK

    K_if   = 25.3403
    K_05if = 18.1115
    n_if   = 9.2453
    V_shift = K_if * (cAMP ** n_if / (K_05if ** n_if + cAMP ** n_if)) - 18.1040
    y_inf = 1.0 / (1.0 + np.exp((Vm + 104.2 - V_shift) / 16.3))
    tau_y = 1.5049 / (np.exp(-(Vm + 590.3) * 0.01094) + np.exp((Vm - 85.1) / 17.2))
    dy = (y_inf - y) / tau_y

    # S.3.3.8 Inward rectifier potassium current
    xk1inf = 1.0 / (1.0 + np.exp(0.070727 * (Vm - EK)))
    IK1 = con.g_K1 * xk1inf * (con.Ko / (con.Ko + 0.228880)) * (Vm - EK)

    # S.3.3.9 Rapidly-activated delayed rectifier potassium current, IKr
    IKr = con.g_Kr * (Vm - EK) * pa * pi
    pa_inf = 1.0 / (1.0 + np.exp(-(Vm + 21.173694) / 9.757086))
    pi_inf = 1.0 / (1.0 + np.exp((Vm + 20.758474 - 4.0) / 19.0))
    tau_pa = 0.699821 / (0.003596 * np.exp(Vm / 15.339290) + 0.000177 * np.exp(-Vm / 25.868423))
    tau_pi = 0.2 + 0.9 / (0.1 * np.exp(Vm / 54.645) + 0.656 * np.exp(Vm / 106.157))
    dpa = (pa_inf - pa) / tau_pa
    dpi = (pi_inf - pi) / tau_pi

    # S.3.3.10 Slowly-activating delayed rectifier potassium current, IKs
    IKs = con.g_Ks * (Vm - EKs) * (xs ** 2)
    xs_inf = 1.0 / (1.0 + np.exp(-(Vm - 20.876040) / 11.852723))
    tau_xs = 1000.0 / (13.097938 / (1.0 + np.exp(-(Vm - 48.910584) / 10.630272)) + np.exp(-Vm / 35.316539))
    dxs = (xs_inf - xs) / tau_xs

    # S.3.3.11 Sodium current, INa
    FNa = ((9.52e-02) * np.exp(-6.3e-2 * (Vm + 34.4)) /(1.0 + 1.66 * np.exp(-0.225 * (Vm + 63.7)))) + 8.69e-2
    hs11 = (1.0 - FNa) * h11 + FNa * j11
    hs15 = (1.0 - FNa) * h15 + FNa * j15

    INa11 = (con.g_Na11 * m11 ** 3 * hs11 * Vm * con.Nao * con.F / (con.E_T * 1000.0) * (np.exp((Vm - ENa) / con.E_T) - 1.0) / (np.exp(Vm / con.E_T) - 1.0))
    
    INa15 = (con.g_Na15 * m15 ** 3 * hs15 * Vm * con.Nao * con.F / (con.E_T * 1000.0) * (np.exp((Vm - con.ENa15) / con.E_T) - 1.0) /(np.exp(Vm / con.E_T) - 1.0))
    
    INa = INa11 + INa15

    m11_inf = 1.0 / (1.0 + np.exp(-(Vm + 36.097331 - 5.0) / 5.0)) ** (1.0 / 3.0)
    h11_inf = 1.0 / (1.0 + np.exp((Vm + 56.0) / 3.0))
    j11_inf = h11_inf
    m15_inf = 1.0 / (1.0 + np.exp(-(Vm + 45.213705) / 7.219547)) ** (1.0 / 3.0)
    h15_inf = 1.0 / (1.0 + np.exp(-(Vm + 62.578120) / (-6.084036)))
    j15_inf = h15_inf

    tau_m11 = 1000.0 * (0.6247e-03 / (0.832 * np.exp(-0.335 * (Vm + 56.7)) + 0.627 * np.exp(0.082 * (Vm + 65.01))) + 0.0000492)
    tau_h11 = 1000.0 * ((3.717e-06 * np.exp(-0.2815 * (Vm + 17.11)) / (1.0 + 0.003732 * np.exp(-0.3426 * (Vm + 37.76)))) + 0.0005977)
    tau_j11 = 1000.0 * ((3.186e-08 * np.exp(-0.6219 * (Vm + 18.8)) / (1.0 + 7.189e-05 * np.exp(-0.6683 * (Vm + 34.07)))) + 0.003556)
    tau_m15 = tau_m11
    tau_h15 = tau_h11
    tau_j15 = tau_j11

    dm15 = (m15_inf - m15) / tau_m15
    dh15 = (h15_inf - h15) / tau_h15
    dj15 = (j15_inf - j15) / tau_j15
    dm11 = (m11_inf - m11) / tau_m11
    dh11 = (h11_inf - h11) / tau_h11
    dj11 = (j11_inf - j11) / tau_j11

    # S.3.3.12 Na+ - K+ pump current, INaK
    INaK = (con.I_NaKmax * con.Ko**1.2 / (con.K_mK**1.2 + con.Ko**1.2)) * \
       (con.Nai**1.3 / (con.K_mNa**1.3 + con.Nai**1.3)) / \
       (1.0 + np.exp(-(Vm - ENa + 120.0) / 30.0))


    # 3.3.13 INaCa (NCX)
    d0 = (1.0 + (con.Cao / con.K_co) * (1.0 + np.exp(con.Q_co * Vm / con.E_T)) + (con.Nao / con.K_1no) * (1.0 + (con.Nao / con.K_2no) * (1.0 + con.Nao / con.K_3no)))
    k43 = con.Nai / (con.K_3ni + con.Nai)
    k41 = np.exp(-con.Q_n * Vm / (2.0 * con.E_T))
    k34 = con.Nao / (con.K_3no + con.Nao)
    k21 = (con.Cao / con.K_co) * np.exp(con.Q_co * Vm / con.E_T) / d0
    k23 = ((con.Nao / con.K_1no) * (con.Nao / con.K_2no) * (1.0 + con.Nao / con.K_3no) * np.exp(-con.Q_n * Vm / (2.0 * con.E_T)) / d0)
    k32 = np.exp(con.Q_n * Vm / (2.0 * con.E_T))
    x1 = k34 * k41 * (k23 + k21) + k21 * k32 * (k43 + k41)
    di = (1.0 + (Casub / con.K_ci) * (1.0 + np.exp(-con.Q_ci * Vm / con.E_T) + con.Nai / con.K_cni) + (con.Nai / con.K_1ni) * (1.0 + (con.Nai / con.K_2ni) * (1.0 + con.Nai / con.K_3ni)))
    k12 = (Casub / con.K_ci) * np.exp(-con.Q_ci * Vm / con.E_T) / di
    k14 = ((con.Nai / con.K_1ni) * (con.Nai / con.K_2ni) * (1.0 + con.Nai / con.K_3ni) * np.exp(con.Q_n * Vm / (2.0 * con.E_T)) / di)
    x2 = k43 * k32 * (k14 + k12) + k41 * k12 * (k34 + k32)
    x3 = k43 * k14 * (k23 + k21) + k12 * k23 * (k43 + k41)
    x4 = k34 * k23 * (k14 + k12) + k21 * k14 * (k34 + k32)

    INaCa = con.K_NaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4)

    # 3.3.14 Sustained inward current , Ist
    I_st = con.g_st * (Vm - con.E_st) * qa * qi
    qa_inf = 1.0 / (1.0 + np.exp(-(Vm + 67.0) / 5.0))
    alpha_qa = 1.0 / (0.15 * np.exp(-Vm / 11.0) + 0.2  * np.exp(-Vm / 700.0))
    beta_qa = 1.0 / (16.0 * np.exp(Vm / 8.0) + 15.0 * np.exp(Vm / 50.0))
    tau_qa = 1.0 / (alpha_qa + beta_qa)
    alpha_qi = 0.15 / (3100.0 * np.exp((Vm + 10.0) / 13.0) + 700.3 * np.exp((Vm + 10.0) / 70.0))
    beta_qi = (0.15 / (95.7 * np.exp(-(Vm + 10.0) / 10.0) + 50.0 * np.exp(-(Vm + 10.0) / 700.0)) + 0.000229 / (1.0 + np.exp(-(Vm + 10.0) / 5.0)))
    qi_inf = alpha_qi / (alpha_qi + beta_qi)
    tau_qi = 1.0 / (alpha_qi + beta_qi)

    dqa = (qa_inf - qa) / tau_qa
    dqi = (qi_inf - qi) / tau_qi

    # 3.3.15 Acetylcholine-activated potassium current,IKACh
    I_KACh = con.C * con.g_KACh_max * (Vm - EK) * w
    # avoid division by zero when CCh = 0
    if con.CCh != 0:
        beta_w = 0.001 * 12.32 / (1.0 + 0.0042 / (con.CCh * 10.0 ** -6))
    else:
        beta_w = 0.0
    alpha_w = 0.001 * 17.0 * np.exp(0.0133 * (Vm + 40.0))
    denom_w = alpha_w + beta_w
    if denom_w != 0:
        w_inf = beta_w / denom_w
        tau_w = 1.0 / denom_w
        a_w = w_inf / tau_w
        b_w = (1.0 - w_inf) / tau_w
    else:
        w_inf = 0.0
        tau_w = np.inf
        a_w = 0.0
        b_w = 0.0
    dw = a_w * (1.0 - w) - b_w * w

    # ============== S.3.4 Ca2+ fluxes in SR ==============
 
    #S.3.4.1 Ryanodine receptor function

    # *************************Add your code here*************************

    # --- WT Parameters ---
    k_om = 0.06   # Transition rate from open to closed state [ms^-1]
    n_Ca = 2.0    # Power of subspace Ca2+ in transition equations
    
    # Equation (1): Ca2+ release flux from the JSR to the subspace
    jSRCarel = con.k_s * OO * (CajSR - Casub)
    
    # Equation (2): SR Ca2+ release termination regulator
    kCaSR = con.MaxSR - (con.MaxSR - con.MinSR) / (1.0 + (con.EC_50SR / CajSR)**con.HSR)
    
    # Calculation of k_oCa based on PKA activity
    k_oCa = con.k_oCa_max * (con.RyR_min - con.RyR_max * PKA * con.n_RyR / (con.k_05Ry * con.n_RyR + PKA ** con.n_RyR) + 1.0)
    
    # Equations (3-4): Transition rates affected by JSR Ca2+ load
    koSRCa = k_oCa / kCaSR
    kiSRCa = con.k_iCa * kCaSR

    # Equations (5-8): RyR 4-state Markov model derivatives
    dR  = (con.k_im * RI - kiSRCa * Casub * R) - (koSRCa * (Casub**n_Ca) * R - k_om * OO)
    dOO = (koSRCa * (Casub**n_Ca) * R - k_om * OO) - (kiSRCa * Casub * OO - con.k_im * S)
    dS  = (kiSRCa * Casub * OO - con.k_im * S) - (k_om * S - koSRCa * (Casub**n_Ca) * RI)
    dRI = (k_om * S - koSRCa * (Casub**n_Ca) * RI) - (con.k_im * RI - kiSRCa * Casub * R)

    # **********************************************************************
    
    
    # S.3.4.2 The rate of calcium uptake (pumping) by the SR, jup
    fPLB = 2.9102 * PLB ** 9.5517 / (0.2763 ** 9.5517 + PLB ** 9.5517) + 0.4998
    j_up = con.P_up * fPLB / (1.0 + con.K_up / Cai)

    # S.3.4.3 Calcium diffusion flux from submembrane space to myoplasm,jCadiff 
    j_Cadiff = (Casub - Cai) / con.tho_difCa

    # S.3.4.4 Calcium flux between network and junctional SR compartments, jtr
    j_tr = (CanSR - CajSR) / con.tho_tr

    # S.3.4.7 Force
    Ve  = 0.0
    dSL = -Ve
    NXB = (SL - con.SL_lo) / 2.0 * 1e3 * (TT + U) * con.Nc
    K_Ca = con.FK0 + con.Fkl * NXB ** con.FN / (con.FK05 ** con.FN + NXB ** con.FN)
    k_l  = con.Fkl / K_Ca
    dA  = con.Fkl * Cai * (1.0 - A - TT - U) - (con.Ff + k_l) * A + (con.Fg0 + con.Fg1 * Ve) * TT
    dTT = con.Ff * A - (con.Fg0 + con.Fg1 * Ve + k_l) * TT + con.Fkl * Cai * U
    dU  = k_l * TT - (con.Fg0 + con.Fg1 * Ve + con.Fkl * Cai) * U

    # S.3.4.5 Natural buffering
    df_TC  = con.Fkl * Cai * (1.0 - A - TT) - k_l * (A + TT)
    df_TMC = con.k_fTMC * Cai * (1.0 - f_TMC - f_TMM) - con.k_bTMC * f_TMC
    df_TMM = con.k_fTMM * con.Mgi * (1.0 - f_TMC - f_TMM) - con.k_bTMM * f_TMM
    df_CMi = con.k_fCM  * Cai   * (1.0 - f_CMi) - con.k_bCM * f_CMi
    df_CMs = con.k_fCM  * Casub * (1.0 - f_CMs) - con.k_bCM * f_CMs

    # *************************Add your code here***************************
    # WT value for calsequestrin association constant
    df_CQ = con.k_fCQ * CajSR * (1.0 - f_CQ) - con.k_bCQ * f_CQ

    # **********************************************************************

    # S.3.4.6 Ca dynamics
    dCai = ((j_Cadiff * con.V_sub - j_up * con.V_nSR) / con.V_i - (con.CM_tot * df_CMi + con.TC_tot * df_TC + con.TMC_tot * df_TMC))
    dCasub = ((-(ICaL + ICaT + Ib_Ca - 2.0 * INaCa) / (2.0 * con.F / 1000.0) + jSRCarel * con.V_jSR) / con.V_sub- j_Cadiff- con.CM_tot * df_CMs)
    
    # *************************Add your code here**************************
   
    # Equation (10): JSR calcium concentration derivative 
    dCajSR = j_tr - jSRCarel - con.CQ_tot * df_CQ
    # **********************************************************************

    dCanSR = j_up - j_tr * con.V_jSR / con.V_nSR

    # ============== 3.5 Membrane potential ==============
    I_tot = (INa + ICaL + ICaT + IKr + IKs + If + It0 + Isus + IK1 + INaK + INaCa + Ib_Na + Ib_K + Ib_Ca + I_KACh)
    dVm = -(I_tot + I_st) / con.C

    # ============== 3.6 Derivatives ==============
    dY = np.zeros(43, dtype=float)
    dY[0]  = dVm
    dY[1]  = dqa
    dY[2]  = dqi
    dY[3]  = ddT
    dY[4]  = dfT
    dY[5]  = dpa
    dY[6]  = dpi
    dY[7]  = dxs
    dY[8]  = dfL12
    dY[9]  = ddL12
    dY[10] = dfL13
    dY[11] = ddL13
    dY[12] = dfCa
    dY[13] = dr
    dY[14] = dq
    dY[15] = dm15
    dY[16] = dh15
    dY[17] = dj15
    dY[18] = dm11
    dY[19] = dh11
    dY[20] = dj11
    dY[21] = dy
    dY[22] = dCai
    dY[23] = dCajSR
    dY[24] = dCanSR
    dY[25] = dCasub
    dY[26] = df_TC
    dY[27] = df_TMC
    dY[28] = df_TMM
    dY[29] = df_CMs
    dY[30] = df_CMi
    dY[31] = df_CQ
    dY[32] = dR
    dY[33] = dOO
    dY[34] = dS
    dY[35] = dRI
    dY[36] = dw
    dY[37] = dcAMP
    dY[38] = dPLB
    dY[39] = dA
    dY[40] = dTT
    dY[41] = dU
    dY[42] = dSL

    return dY
