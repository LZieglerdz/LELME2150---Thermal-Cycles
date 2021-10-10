"""
LELME2150 - Thermal cycles
Homework 2 - Exergy and combustion

Signature of the final gas turbine function

@author: Antoine Laterre
@date: September 22, 2021
"""

#
#===IMPORT PACKAGES============================================================
#

import CoolProp.CoolProp as CP
import numpy as np
# import re      #try to use regex to get composition of fuel from chemical formula: C4H10 -> C:4, H:10;

#
#===GLOBAL VARIABLES===========================================================
#
T_S = 273.153         # [K]
fuel = 'CH4'
LHV = 50150e3         # J/kg_CH4
# # air composition -- non trivial
# N2_conc = .7667       # %mol        trivial value: .79
# O2_conc = .2064       # %mol        trivial value: .21
# H2O_conc = .0170      # %mol
# Ar_conc = .0095       # %mol
# CO2_conc = .0004      # %mol
# air composition -- trivial
N2_conc = .79         # %mol        trivial value: .79
O2_conc = .21         # %mol        trivial value: .21
H2O_conc = 0          # %mol
Ar_conc = 0           # %mol
CO2_conc = 0          # %mol

Mm_N2 = CP.PropsSI("MOLARMASS", "N2")        # kg/mol
Mm_O2 = CP.PropsSI("MOLARMASS", "O2")        # kg/mol
Mm_H2O = CP.PropsSI("MOLARMASS", "H2O")        # kg/mol
Mm_H2 = CP.PropsSI("MOLARMASS", "H2")        # kg/mol
Mm_Ar = CP.PropsSI("MOLARMASS", "Ar")        # kg/mol
Mm_CO2 = CP.PropsSI("MOLARMASS", "CO2")        # kg/mol
Mm_fuel = CP.PropsSI("MOLARMASS", fuel)        # kg/mol
Mm_C = Mm_CO2-Mm_O2        # kg/mol
Mm_H = .5*Mm_H2         # kg/mol
Mm_O = .5*Mm_O2      # kg/mol

comp = ["N2", "O2", "H2O", "Ar", "CO2"]
conc = [N2_conc, O2_conc, H2O_conc,Ar_conc, CO2_conc]
Mm   = [Mm_N2, Mm_O2, Mm_H2O, Mm_Ar, Mm_CO2]
Mm_air = np.dot(Mm, conc)
#
#===STOECHIOMETRY===========================================================
#
def get_stoech(fuel, z, y, x, T_2, T_3, p_2, p_3, R_Star, iter, w):
    #  Combustion supposée complète
    #  C_zH_yO_x + w(O_2 + 3.76 N_2) --> a_0 O_2 + a_2 CO_2 + b_2 H_2O + 3.76w N_2
    # comb stoechiometric -> w = z + (y-2*x)/4
    # a_2 = z
    # 2*b_2 = y
    # 2*a_0 + 2*a_2 + b_2 = x + 2*w

    # ma1 = ( (Mm_O2 + 3.76*Mm_N2) * (z+ (y-2*x)/4) ) / (z*Mm_C + y*Mm_H + x*Mm_O)

    a_0 = (x + 2*w - (y/2) - 2*z)/2

    sum_R = getMeanCp(p_2, p_3, T_2, T_3, R_Star, fuel) + w*getMeanCp(p_2, p_3, T_2, T_3, R_Star, 'air')
    sum_P = a_0*getMeanCp(p_2, p_3, T_2, T_3, R_Star,'O2') + z*getMeanCp(p_2, p_3, T_2, T_3, R_Star, 'CO2') + (y/2)*getMeanCp(p_2, p_3, T_2, T_3, R_Star, 'H2O') + 3.76*w*getMeanCp(p_2, p_3, T_2, T_3, R_Star, 'N2')
    delta_T = (LHV + (T_2-T_S)*sum_R)/sum_P
    error = delta_T - (T_3-T_S)
    print(iter, w, error)
    if iter == 0:
        print('The function get_stoech did not converge.')
        return(w)
    if np.abs(error) <= 1e-6:
        return(w)
    else:
        return(get_steoch(fuel, z, y, x, T_2, T_3, p_2, p_3, R_Star,iter-1, w+.002))

def get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3, h_2, R_Star, iter, lam_est):
    #  Combustion supposée complète
    #  C_zH_yO_x + w(O_2 + 3.76 N_2) --> a_0 O_2 + a_2 CO_2 + b_2 H_2O + 3.76w N_2
    # comb stoechiometric -> w = z + (y-2*x)/4
    # a_2 = z
    # 2*b_2 = y
    # 2*a_0 + 2*a_2 + b_2 = x + 2*w
    w = lam_est*(z + (y-2*x)/4)
    a_0 = (x + 2*w - (y/2) - 2*z)/2
    m_fuel = 1*Mm_fuel
    m_air = w*Mm_air
    m_flue = a_0*Mm_O2 + z*Mm_CO2 + y*Mm_H2O/2 + 3.76*w*Mm_N2
    mol_tot = (3.76*w + a_0 + y/2 + z)
    conc_flue = [3.76*w/mol_tot, a_0/mol_tot, y/2/mol_tot, 0, z/mol_tot]

    ma1 = ( (Mm_O2 + 3.76*Mm_N2) * (z+ (y-2*x)/4) ) / (z*Mm_C + y*Mm_H + x*Mm_O)

    h_3 = h_2 + (m_fuel*getMeanCp(p_2, p_3, T_2, T_3, R_Star,fuel) + m_air*getMeanCp(p_2, p_3, T_2, T_3, R_Star,'air') + m_flue*getMeanCp2(p_2, p_3, T_2, T_3, R_Star, comp,conc_flue))*(T_3-T_2)#/ (np.dot(Mm, conc_flue))
    lam = (LHV-h_3)/ (ma1*(h_3-h_2))
    # print(iter, lam)
    if iter == 0:
        print('The function get_lambda did not converge.')
        return(-1)
    if np.abs(lam-lam_est) <= 1e-6 :
        return(ma1, lam)
    else:
        return( get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3, h_2, R_Star, iter-1, lam) )

def get_lambda2(fuel, z, y, x, T_2, T_3, p_2, p_3, h_2, R_Star, iter, lam_est):
    #  Combustion supposée complète
    #  C_zH_yO_x + w(O_2 + 3.76 N_2) --> a_0 O_2 + a_2 CO_2 + b_2 H_2O + 3.76w N_2
    # comb stoechiometric -> w = z + (y-2*x)/4
    # a_2 = z
    # 2*b_2 = y
    # 2*a_0 + 2*a_2 + b_2 = x + 2*w
    w = lam_est*(z + (y-2*x)/4)
    a_0 = (x + 2*w - (y/2) - 2*z)/2
    m_fuel = 1*Mm_fuel
    m_air = w*Mm_air
    m_flue = a_0*Mm_O2 + z*Mm_CO2 + y*Mm_H2O/2 + 3.76*w*Mm_N2
    mol_tot = (3.76*w + a_0 + y/2 + z)
    conc_flue = [3.76*w/mol_tot, a_0/mol_tot, y/2/mol_tot, 0, z/mol_tot]

    ma1 = ( (Mm_O2 + 3.76*Mm_N2) * (z+ (y-2*x)/4) ) / (z*Mm_C + y*Mm_H + x*Mm_O)
    sum_Cp = (m_fuel*getMeanCp(p_2, p_3, T_2, T_3, R_Star,fuel) + m_air*getMeanCp(p_2, p_3, T_2, T_3, R_Star,'air') + m_flue*getMeanCp2(p_2, p_3, T_2, T_3, R_Star, comp,conc_flue))
    h_3 = sum_Cp*(T_3-T_S)#/(np.dot(Mm, conc_flue))
    lam = (LHV-h_3)/ (ma1*sum_Cp*(T_3-T_2))
    # print(iter, lam)
    if iter == 0:
        print('The function get_lambda did not converge.')
        return(-1)
    if np.abs(lam-lam_est) <= 1e-6 :
        return(ma1, lam)
    else:
        return( get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3, h_2, R_Star, iter-1, lam) )

#
#===MEAN CP====================================================================
#

def getCpMix(T, p, mix_comp, mix_conc):
    sum = 0
    for i in range(len(mix_comp)):
         sum += (mix_conc[i]*CP.PropsSI('CPMASS', 'T', T, 'P', p, mix_comp[i]) )
    return(sum)

def getMeanCp(p_in, p_out, T_in, T_out, R_Star, str):               #renvoie le cp massique moyen
    T_min = np.min([T_in,T_out])
    T_max = np.max([T_in,T_out])
    p_min = np.min([p_in,p_out])
    p_max = np.max([p_in,p_out])
    if str == 'air':
        mix_comp = comp
        mix_conc = conc
        cp = getCpMix(T_in, p_in, mix_comp, mix_conc)
    elif str == 'flue':
        mix_comp = comp                                     # a modif
        mix_conc = conc                                     # a modif
        cp = getCpMix(T_in, p_in, mix_comp, mix_conc)
    else:
        mix_comp = [str]
        mix_conc = [1]
        cp = getCpMix(T_in, p_in, mix_comp, mix_conc)
    if T_min == T_max:
        return(cp)
    n=100
    rangeT = np.linspace(T_min,T_max,n)
    rangeP = np.ones(n)*p_min
    p = p_min
    for i in np.arange(1,n):
        gamma = cp/(cp - i*R_Star)
        p = p_min*(rangeT[i]/T_min)**((gamma-1)/gamma)
        cp += getCpMix(rangeT[i], rangeP[i], mix_comp, mix_conc)
    return(cp/n)

def getMeanCp2(p_in, p_out, T_in, T_out, R_Star, mix_comp, mix_conc):               #renvoie le cp massique moyen
    T_min = np.min([T_in,T_out])
    T_max = np.max([T_in,T_out])
    p_min = np.min([p_in,p_out])
    p_max = np.max([p_in,p_out])
    cp = getCpMix(T_in, p_in, mix_comp, mix_conc)

    if T_min == T_max:
        return(cp)
    n=100
    rangeT = np.linspace(T_min,T_max,n)
    rangeP = np.ones(n)*p_min
    p = p_min
    for i in np.arange(1,n):
        gamma = cp/(cp - i*R_Star)
        p = p_min*(rangeT[i]/T_min)**((gamma-1)/gamma)
        cp += getCpMix(rangeT[i], rangeP[i], mix_comp, mix_conc)
    return(cp/n)


#
#===POLYTROPIC TEMPERATURE=====================================================
#
def getPolytropicTemp(p_in, p_out, T_in, T_out, R_Star, eta_pi, iter, str):
    if iter < 0:
        print("Function does not converge.")
        return(-1)
    # T = T_in*(p_out/p_in)**(R_Star/(eta_pi*getMeanCp(p_in,p_out,T_in,T_out, R_Star, str )))
    T = T_in*(p_out/p_in)**(R_Star/(eta_pi*getMeanCp(p_in,p_out,T_in,T_out, R_Star, str )))
    if np.abs(T_out-T) < 1e-12:
        return(T)
    else:
        return(getPolytropicTemp(p_in, p_out, T_in, T, R_Star, eta_pi, iter-1, str))


#
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

# IN:
#   - P_e [W]: net power output
#   - options: tuple containing the parametric values of the GT
#       o p_1 [Pa], T_1 [K]: inlet pressure and temperature (ambient)
#       o r [-], eta_pi_c [-]: compression ratio, compressor polytropic efficiency
#       o T_3 [K], k_cc [-]: combustor outlet temperature and pressure losses coeffcient
#       o eta_pi_t [-]: turbine polytropic efficiency
#       o k_mec [-]: shaft losses
#   - display: bool to choose to plot the T-s & h-s diagrams and the energy and exergy pie charts (True or False)
# OUT: tuple containing
#   - DAT: tuple containing the cycle state data
#      o p [Pa]: tuple containing the pressure at each state
#      o T [T]: tuple containing the temperature at each state
#      o h [J/kg]: tuple containing the enthalpy at each state
#      o s [J/kg/K]: tuple containing the entropy at each state
#      o e [J/kg]: tuple containing the exergy at each state
#   - COMBUSTION: tuple containing the combustion parameters
#      o LHV [J/kg]: the fuel Lower Heating Value
#      o e_c [J/kg]: the fuel exergy
#      o excess_air [-]: the excess air of the combustion
#      o cp_gas [J/kg/K]: the flue gas specific heat capacity at the combustor outlet
#      o gas_prop: list containing the proportion of ['CO2','H2O','N2','O2'] in the flue gas respectively
#   - MASSFLOW: tuple containing the massflow rates
#      o dot_m_a [kg/s]: mass flow rate of air
#      o dot_m_f [kg/s]: mass flow rate of fuel
#      o dot_m_g [kg/s]: mass flow rate of flue gas
#   - ETA: tuple containing the efficiencies (see text book pp. 113-125)
#       o eta_cyclen [-]: cycle energy efficiency
#       o eta_toten [-]: overall energy efficiency
#       o eta_cyclex [-]: cycle exergy efficiency
#       o eta_totex [-]: overall exergy efficiency
#       o eta_rotex [-]: compressor-turbine exergy efficiency
#       o eta_combex [-]: combustion exergy efficiency
#   - DATEN: tuple containing the energy losses
#       o loss_mec [W]: mechanical energy losses
#       o loss_ech [W]: exhaust energy losses
#   - DATEX: tuple containing the exergy losses
#       o loss_mec [W]: mechanical energy losses
#       o loss_rotex [W]: compressor-turbine exergy losses
#       o loss_combex [W]: combustion exergy losses
#       o loss_echex [W]: exhaust exergy losses
#   - FIG: tuple containing the figures to be diplayed ((0,0,0,0) if display = False)
#       o fig_pie_en: pie chart of energy losses
#       o fig_pie_ex: pie chart of exergy losses
#       o fig_Ts_diagram: T-s diagram of the GT cycle
#       o fig_hs_diagram: h-s diagram of the GT cycle
# REMARK:
#   - You can adapt the numerical application of the text book
#   in page pp. 125 to check your results.
# - As improvements, other options can be added (not expected before HW4):
#       o air [n.a.], air_prop [-]: air composition (types: list of str, list of float)
#       o NTU [-]: number of transfer unit in case of a recuperator
#       o fuel [n.a.]: type of fuel  (type: str), default is 'methane'
def gas_turbine(P_e,options,display):
    # Process input variables--------------------------------------------------
    p_1,T_1,r,eta_pi_c,T_3,k_cc,eta_pi_t,k_mec = options

    # Replace with your model--------------------------------------------------
    p_1, p_2, p_3, p_4 = p_1, 0, 0, 0
    T_1, T_2, T_3, T_4 = T_1, 0, T_3, 0
    h_1, h_2, h_3, h_4 = 0, 0, 0, 0
    s_1, s_2, s_3, s_4 = 0, 0, 0, 0
    e_1, e_2, e_3, e_4 = 0, 0, 0, 0
    LHV,e_c,excess_air,cp_gas,gas_prop = 0, 0, 0, 0, 0
    dotm_a,dotm_f,dotm_g = 0, 0, 0
    eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex = 0, 0, 0, 0, 0, 0
    loss_mec,loss_ech = 0, 0
    loss_mec,loss_rotex,loss_combex,loss_echex = 0, 0, 0, 0
    fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram = 0, 0, 0, 0

    R_Star = 0
    for i in range(len(comp)):
        R_Star += conc[i]*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, comp[i])/Mm[i] # [J/kgK]


    #set references state
    for i in comp:
        Dmolar = CP.PropsSI("Dmolar", "T", T_S, "P", p_1, i)
        CP.set_reference_state(i, T_S, Dmolar, 0, 0)


    # State 1 -- 4->1: Isobar Heat Rejection
    s_1 = .79*CP.PropsSI('S','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('S','T', T_1,'P', p_1, "O2")
    h_1 = .79*CP.PropsSI('H','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('H','T', T_1,'P', p_1, "O2")
    e_1 = h_1 - T_1*s_1

    # State 2 -- 1->2: Polytropic Compression
    p_2 = r*p_1
    T_2 = getPolytropicTemp(p_1, p_2, T_1, T_1, R_Star, eta_pi_c, 100, 'air')
    h_2 = getMeanCp(p_1,p_2,T_1,T_2, R_Star,'air')*(T_2-T_1)
    s_2 = s_1 + getMeanCp(p_1,p_2,T_1,T_2, R_Star, 'air')*np.log(T_2/T_1) - R_Star*np.log(p_2/p_1)
    e_2 = h_2 - T_2*s_2

    # State 3 -- 2->3: Isobar Combustion
    p_3 = k_cc*p_2
    s_3 = .79*CP.PropsSI('S','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('S','T', T_3,'P', p_3, "O2")
    h_3 = .79*CP.PropsSI('H','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('H','T', T_3,'P', p_3, "O2")
    e_3 = h_3 - T_3*s_3

    # State 4 -- 3->4: Polytropic Expansion
    p_4 = p_1
    T_4 = getPolytropicTemp(p_3, p_4, T_3, T_3, R_Star, 1/eta_pi_t, 100, 'flue')                      # CHANGE TRUE TO FALSE WHEN CPGAS IMPLEMENTED
    h_4 = h_3 + getMeanCp(p_3,p_4,T_4,T_3, R_Star, 'flue')*(T_4-T_3)                                  # CHANGE TRUE TO FALSE WHEN CPGAS IMPLEMENTED
    s_4 = s_3 + getMeanCp(p_3,p_4,T_3,T_4, R_Star, 'flue')*np.log(T_4/T_3) - R_Star*np.log(p_4/p_3)   # CHANGE TRUE TO FALSE WHEN CPGAS IMPLEMENTED
    e_4 = h_4 - T_4*s_4

    # # Efficiency:
    # W = (h_3-h_4)*eta_mec_t - (h_2-h_1)/eta_mec_c
    # Q = h_3-h_2
    # eta_cyclen = W/Q
    # eta_mec = 1
    # eta_gen = 1
    # eta_toten = eta_cyclen*eta_mec*eta_gen

    # Process output variables - do not modify---------------------------------
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    s = (s_1, s_2, s_3, s_4)
    h = (h_1, h_2, h_3, h_4)
    e = (e_1, e_2, e_3, e_4)
    DAT = (p,T,h,s,e)
    COMBUSTION = (LHV,e_c,excess_air,cp_gas,gas_prop)
    MASSFLOW = (dotm_a,dotm_f,dotm_g)
    ETA = (eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex)
    DATEN = (loss_mec,loss_ech)
    DATEX = (loss_mec,loss_rotex,loss_combex,loss_echex)
    FIG = (fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram)
    out = (ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG)
    return out
