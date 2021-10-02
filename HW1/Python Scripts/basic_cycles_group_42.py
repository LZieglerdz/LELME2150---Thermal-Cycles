"""
LELME2150 - Thermal cycles
Homework 1 - Basic cycles

Signatures of the functions for the homework

@author: Antoine Laterre
@date: September 10, 2021
"""

#
#===IMPORT PACKAGES============================================================
#

import CoolProp.CoolProp as CP
import numpy as np

#
#=== MEAN CP ==========================================
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
#=== BRAYTON CYCLE - GAS TURBINE ==========================================
#
def gas_turbine(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t):

    p_1, p_2, p_3, p_4 = p_1, p_2, p_3, p_1     # [Pa]
    T_1, T_2, T_3, T_4 = T_1, 0, T_3, 0         # [K]
    s_1, s_2, s_3, s_4 = 0, 0, 0, 0             # [J/kgK]
    h_1, h_2, h_3, h_4 = 0, 0, 0, 0             # [J/kg]
    eta_en = 0

    R_Star = .79*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, "N2")/CP.PropsSI('MOLAR_MASS','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, "O2")/CP.PropsSI('MOLAR_MASS','T', T_1,'P', p_1, "O2") # [J/kgK]

    #set references state
    DmolarN2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "N2")
    CP.set_reference_state("N2", T_1, DmolarN2, 0, 0)
    DmolarO2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "O2")
    CP.set_reference_state("O2", T_1, DmolarO2, 0, 0)

    # State 1 -- 4->1: Isobar Heat Rejection
    s_1 = .79*CP.PropsSI('S','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('S','T', T_1,'P', p_1, "O2")
    h_1 = .79*CP.PropsSI('H','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('H','T', T_1,'P', p_1, "O2")

    # State 2 -- 1->2: Polytropic Compression
    T_2 = getPolytropicTemp(p_1, p_2, T_1, T_1, R_Star, eta_pi, 100)
    h_2 = getMeanCp(p_1,p_2,T_1,T_2, R_Star)*(T_2-T_1)
    s_2 = s_1 + getMeanCp(p_1,p_2,T_1,T_2, R_Star)*np.log(T_2/T_1) - R_Star*np.log(p_2/p_1)

    # State 3 -- 2->3: Isobar Combustion
    s_3 = .79*CP.PropsSI('S','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('S','T', T_3,'P', p_3, "O2")
    h_3 = .79*CP.PropsSI('H','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('H','T', T_3,'P', p_3, "O2")

    # State 4 -- 3->4: Polytropic Expansion
    T_4 = getPolytropicTemp(p_3, p_4, T_3, T_3, R_Star, 1/eta_pi, 100)
    h_4 = h_3 + getMeanCp(p_3,p_4,T_4,T_3, R_Star)*(T_4-T_3)
    s_4 = s_3 + getMeanCp(p_3,p_4,T_3,T_4, R_Star)*np.log(T_4/T_3) - R_Star*np.log(p_4/p_3)

    # Efficiency
    W = (h_3-h_4)*eta_mec_t - (h_2-h_1)/eta_mec_c
    Q = h_3-h_2
    eta_cyclen = W/Q
    eta_mec = 1
    eta_gen = 1
    eta_en = eta_cyclen*eta_mec*eta_gen


    # Final outputs - do not modify
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    s = (s_1, s_2, s_3, s_4)
    h = (h_1, h_2, h_3, h_4)
    out = (p,T,s,h,eta_en)
    # print('     | p\t| T \t\t| s-s1\t\t| h-h1')
    # print('     | [kPa]\t| [K] \t\t| [kJ/kg.K]\t| [kJ/kg]')
    # print(75*'-')
    # print('  1  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p[0]/1000, T[0], s[0]/1000, h[0]/1000))
    # print('  2  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p[1]/1000, T[1], s[1]/1000, h[1]/1000))
    # print('  3  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p[2]/1000, T[2], s[2]/1000, h[2]/1000))
    # print('  4  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p[3]/1000, T[3], s[3]/1000, h[3]/1000))
    # print(75*'-')
    # print('eta_en : %.3f' %eta_en)
    # print('eta_pi : %.2f' %eta_pi)
    # print('eta_mec_c : %.2f' %eta_mec_c)
    # print('W : %.2f [kj/kg]' %(W/1000))
    # print('Q : %.2f [kj/kg]' %(Q/1000))
    # print(75*'_')
    # print('\n')
    return out

#
#=== RANKINE CYCLE - STEAM TURBINE ==========================================
#

def steam_turbine(T_1,p_3,T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump):
   #State 1: Saturated liquid
    T_1 = T_1
    x_1 = 0
    p_1 = CP.PropsSI('P','T',T_1,'Q',0,'Water')
    h_1 = CP.PropsSI('H','P',p_1,'Q',0,'Water')
    s_1 = CP.PropsSI('S','P',p_1,'Q',0,'Water')

    #State 2: Sub-cooled water (1->2 supposed isothermal, 2->3 supposed isobaric)
    p_2 = p_3
    rho = CP.PropsSI('D','T',T_1,'P',p_1,'Water')
    T_2 = T_1 + (p_2-p_1)*(1/eta_pump -1)/(rho*CP.PropsSI('C','P',p_1,'T',T_1,'Water'))
    h_2 = CP.PropsSI('H','P',p_2,'T',T_2,'Water')
    s_2 = CP.PropsSI('S','P',p_2,'T',T_2,'Water')
    x_2 = CP.PropsSI('Q','P',p_2,'T',T_2,'Water')

    #State 3:
    p_3 = p_3
    T_3 = T_3
    h_3 = CP.PropsSI('H','P',p_3,'T',T_3,'Water')
    s_3 = CP.PropsSI('S','P',p_3,'T',T_3,'Water')
    x_3 = CP.PropsSI('Q','P',p_3,'T',T_3,'Water')

    #State 4_si: (Isentropic expansion from 3->4_si)
    p_4_si = p_1 # 4_si->1 : isobaric condensation
    T_4_si = T_1
    s_4_si = s_3
    s_4_si_satLiq = CP.PropsSI('S','P',p_4_si,'Q',0,'Water')
    s_4_si_satVap = CP.PropsSI('S','P',p_4_si,'Q',1,'Water')
    x_4_si = (s_4_si - s_4_si_satLiq)/(s_4_si_satVap - s_4_si_satLiq)
    h_4_si_satLiq = CP.PropsSI('H','P',p_4_si,'Q',0,'Water')
    h_4_si_satVap = CP.PropsSI('H','P',p_4_si,'Q',1,'Water')
    h_4_si = x_4_si*h_4_si_satVap + (1-x_4_si)*h_4_si_satLiq

    #State 4: (Adiabatic expansion from 3->4)
    p_4 = p_1
    T_4 = T_1
    h_4 = h_3 + eta_is_t*(h_4_si - h_3)
    h_4_satLiq = CP.PropsSI('H','P',p_4,'Q',0,'Water')
    h_4_satVap = CP.PropsSI('H','P',p_4,'Q',1,'Water')
    x_4 = (h_4 - h_4_satLiq)/(h_4_satVap - h_4_satLiq)
    s_4_satLiq = CP.PropsSI('S','P',p_4,'Q',0,'Water')
    s_4_satVap = CP.PropsSI('S','P',p_4,'Q',1,'Water')
    s_4 = x_4*s_4_satVap + (1-x_4)*s_4_satLiq

    #Mechanical work of the turbine:
    wm_t = (h_4 - h_3)

    #Overall efficiency of the cycle:
    eta_th = (h_3 - h_4)/(h_3 - h_2)

    eta_en = eta_mec_t*eta_th*eta_gen

    #Water flow in the cycle:
    dot_m_v = P_el/(eta_mec_t*(h_3-h_4-h_2+h_1))

    #Mass flow of coal used:
    dot_m_c = dot_m_v*(h_3 - h_2)/(eta_gen*LHV)

    #Minimum water flow used at the condenser:
    dT_w = 8 # inlet temperature of 8°C when entering the condenser
    dot_m_w = dot_m_v*(h_4 - h_1)/(CP.PropsSI('C','P',101325,'T',288.15,'Water')*dT_w)

    ##?
    dot_m_f = dot_m_v


    # Final outputs - do not modify
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    s = (s_1, s_2, s_3, s_4)
    h = (h_1, h_2, h_3, h_4)
    x = (x_1, x_2, x_3, x_4)
    out = (p,T,s,h,x,eta_en,dot_m_f)
    return out
