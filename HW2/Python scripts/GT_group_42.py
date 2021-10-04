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

#
#===MEAN CP====================================================================
#
def getMeanCp(p_in, p_out, T_in, T_out, R_Star):
    T_min = np.min([T_in,T_out])
    T_max = np.max([T_in,T_out])
    p_min = np.min([p_in,p_out])
    p_max = np.max([p_in,p_out])
    cp = .79*CP.PropsSI('CPMASS','T', T_min,'P', p_min, "N2") + .21*CP.PropsSI('CPMASS','T', T_min,'P', p_min, "O2")
    if T_min == T_max:
        return(cp)
    n=100
    rangeT = np.linspace(T_min,T_max,n)
    rangeP = np.ones(n)*p_min
    p = p_min
    for i in np.arange(1,n):
        gamma = cp/(cp - i*R_Star)
        p = p_min*(rangeT[i]/T_min)**((gamma-1)/gamma)
        cp += .79*CP.PropsSI('CPMASS','T', rangeT[i],'P', rangeP[i], "N2") + .21*CP.PropsSI('CPMASS','T', rangeT[i],'P', p, "O2")
    return(cp/n)
#
#===POLYTROPIC TEMPERATURE=====================================================
#
def getPolytropicTemp(p_in, p_out, T_in, T_out, R_Star, eta_pi, iter):
    if iter < 0:
        print("Function does not converge.")
        return(-1)
    T = T_in*(p_out/p_in)**(R_Star/(eta_pi*getMeanCp(p_in,p_out,T_in,T_out, R_Star)))
    if np.abs(T_out-T) < 1e-12:
        return(T)
    else:
        return(getPolytropicTemp(p_in, p_out, T_in, T, R_Star, eta_pi, iter-1))


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

    R_Star = .79*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, "N2")/CP.PropsSI('MOLAR_MASS','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, "O2")/CP.PropsSI('MOLAR_MASS','T', T_1,'P', p_1, "O2") # [J/kgK]

    #set references state
    DmolarN2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "N2")
    CP.set_reference_state("N2", T_1, DmolarN2, 0, 0)
    DmolarO2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "O2")
    CP.set_reference_state("O2", T_1, DmolarO2, 0, 0)

    # State 1 -- 4->1: Isobar Heat Rejection
    s_1 = 0 #.79*CP.PropsSI('S','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('S','T', T_1,'P', p_1, "O2")    #Les deux lignes s_1 et h_1 sont useless car on
    h_1 = 0 #.79*CP.PropsSI('H','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('H','T', T_1,'P', p_1, "O2")    #considère déjà que ce sont des points de référence
    e_1 = h_1 - T_1*s_1                                                                                   #et qu'ils sont donc nuls.

    # State 2 -- 1->2: Polytropic Compression
    p_2 = r*p_1
    T_2 = getPolytropicTemp(p_1, p_2, T_1, T_1, R_Star, eta_pi_c, 100)
    h_2 = getMeanCp(p_1,p_2,T_1,T_2, R_Star)*(T_2-T_1)
    s_2 = s_1 + getMeanCp(p_1,p_2,T_1,T_2, R_Star)*np.log(T_2/T_1) - R_Star*np.log(p_2/p_1)
    e_2 = h_2 - T_2*s_2

    # State 3 -- 2->3: Isobar Combustion
    p_3 = p_2
    s_3 = .79*CP.PropsSI('S','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('S','T', T_3,'P', p_3, "O2")
    h_3 = .79*CP.PropsSI('H','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('H','T', T_3,'P', p_3, "O2")
    e_3 = h_3 - T_3*s_3

    # State 4 -- 3->4: Polytropic Expansion
    p_4 = p_3/r
    T_4 = getPolytropicTemp(p_3, p_4, T_3, T_3, R_Star, 1/eta_pi_t, 100)
    h_4 = h_3 + getMeanCp(p_3,p_4,T_4,T_3, R_Star)*(T_4-T_3)
    s_4 = s_3 + getMeanCp(p_3,p_4,T_3,T_4, R_Star)*np.log(T_4/T_3) - R_Star*np.log(p_4/p_3)
    e_4 = h_4 - T_4*s_4

    # # Efficiency
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
