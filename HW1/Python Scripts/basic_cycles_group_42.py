"""
LELME2150 - Thermal cycles
Homework 1 - Basic cycles

Signatures of the functions for the homework

@author: Ziegler de Ziegleck aùf Rheingrüb Laurent, Boterberg Dimitri
@date: September 22, 2021
"""

#
#===IMPORT PACKAGES============================================================
#

import CoolProp.CoolProp as CP
import numpy as np


#
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
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
    cp12 = .79*CP.PropsSI('CPMASS','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('CPMASS','T', T_1,'P', p_1, "O2")
    T_2 = T_1*(p_2/p_1)**(R_Star/(eta_pi*cp12))

    delta_T = 10
    k=20
    n=100
    while k>0 and delta_T > 1e-3 :
        k -= 1
        T_Old = T_2
        range = np.linspace(T_1,T_2,n)
        for i in range:
            cp12 += .79*CP.PropsSI('CPMASS','T', i,'P', p_2, "N2") + .21*CP.PropsSI('CPMASS','T', i,'P', p_2, "O2")
        cp12 = cp12/n
        T_2 = T_1*(p_2/p_1)**(R_Star/(eta_pi*cp12))
        delta_T = np.abs(T_2-T_Old)

    h_2 = .79*CP.PropsSI('H','T', T_2,'P', p_2, "N2") + .21*CP.PropsSI('H','T', T_2,'P', p_2, "O2")
    s_2 = .79*CP.PropsSI('S','T', T_2,'P', p_2, "N2") + .21*CP.PropsSI('S','T', T_2,'P', p_2, "O2")

    # State 3 -- 2->3: Isobar Combustion
    s_3 = .79*CP.PropsSI('S','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('S','T', T_3,'P', p_3, "O2")
    h_3 = .79*CP.PropsSI('H','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('H','T', T_3,'P', p_3, "O2")

    # State 4 -- 3->4: Polytropic Expansion
    cp34 = .79*CP.PropsSI('CPMASS','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('CPMASS','T', T_3,'P', p_3, "O2")
    T_4 = T_3*(p_4/p_3)**(R_Star/(eta_pi*cp34))

    delta_T = 10
    k=20
    n=100
    while k>0 and delta_T > 1e-3 :
        k -= 1
        T_Old = T_4
        range = np.linspace(T_3,T_4,n)
        for i in range:
            cp34 += .79*CP.PropsSI('CPMASS','T', i,'P', p_4, "N2") + .21*CP.PropsSI('CPMASS','T', i,'P', p_4, "O2")
        cp34 = cp34/n
        T_4 = T_3*(p_4/p_3)**(R_Star*eta_pi/(cp34))
        delta_T = np.abs(T_4-T_Old)


    h_4 = .79*CP.PropsSI('H','T', T_4,'P', p_4, "N2") + .21*CP.PropsSI('H','T', T_4,'P', p_4, "O2")
    s_4 = .79*CP.PropsSI('S','T', T_4,'P', p_4, "N2") + .21*CP.PropsSI('S','T', T_4,'P', p_4, "O2")

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
    return out

#
#===RANKINE CYCLE - TO BE IMPLEMENTED==========================================
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
    wm_t = (h_3 - h_4)

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

    ##
    dot_m_f = dot_m_v

    # Final outputs - do not modify
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    s = (s_1, s_2, s_3, s_4)
    h = (h_1, h_2, h_3, h_4)
    x = (x_1, x_2, x_3, x_4)
    out = (p,T,s,h,x,eta_en,dot_m_f)
    return out
