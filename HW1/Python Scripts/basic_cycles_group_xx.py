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
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

# def gas_turbine(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t):
#
#
#     p_1, p_2, p_3, p_4 = p_1, p_2, p_3, p_1  # [Pa]
#     T_1, T_2, T_3, T_4 = T_1, 0, T_3, 0    # [K]
#     s_1, s_2, s_3, s_4 = 0, 0, 0, 0        # [J/kgK]
#     h_1, h_2, h_3, h_4 = 0, 0, 0, 0        # [J/kg]
#
#     eta_en = 0
#
#     cp = 1.005e3          # [J/kgK]
#
#     #set references state
#     DmolarN2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "N2")
#     CP.set_reference_state("N2", T_1, DmolarN2, 0, 0)
#     DmolarO2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "O2")
#     CP.set_reference_state("O2", T_1, DmolarO2, 0, 0)
#
#     # State 1
#     s_1 = .79*CP.PropsSI('S','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('S','T', T_1,'P', p_1, "O2")
#     h_1 = .79*CP.PropsSI('H','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('H','T', T_1,'P', p_1, "O2")
#     # State 2
#     # s_2 = s_1
#     T_2 = T_1*(p_2/p_1)**(.4/(1.4*eta_pi))
#     h_2 = .79*CP.PropsSI('H','T', T_2,'P', p_2, "N2") + .21*CP.PropsSI('H','T', T_2,'P', p_2, "O2")
#     s_2 = .79*CP.PropsSI('S','T', T_2,'P', p_2, "N2") + .21*CP.PropsSI('S','T', T_2,'P', p_2, "O2")
#     # T_2 = .79*CP.PropsSI('T','S', s_2,'P', p_2, "N2") + .21*CP.PropsSI('T','S', s_2,'P', p_2, "O2")
#     h_2 = .79*CP.PropsSI('H','S', s_2,'P', p_2, "N2") + .21*CP.PropsSI('H','S', s_2,'P', p_2, "O2")
#     # State 3
#     s_3 = .79*CP.PropsSI('S','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('S','T', T_3,'P', p_3, "O2")
#     h_3 = .79*CP.PropsSI('H','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('H','T', T_3,'P', p_3, "O2")
#     # State 4
#     # s_4 = s_3
#     T_4 = T_3*(p_4/p_3)**(.4/(1.4*eta_pi))
#     h_4 = .79*CP.PropsSI('H','T', T_4,'P', p_4, "N2") + .21*CP.PropsSI('H','T', T_4,'P', p_4, "O2")
#     s_4 = .79*CP.PropsSI('S','T', T_4,'P', p_4, "N2") + .21*CP.PropsSI('S','T', T_4,'P', p_4, "O2")
#     # T_4 = .79*CP.PropsSI('T','S', s_4,'P', p_4, "N2") + .21*CP.PropsSI('T','S', s_4,'P', p_4, "O2")
#     # h_4 = .79*CP.PropsSI('H','S', s_4,'P', p_4, "N2") + .21*CP.PropsSI('H','S', s_4,'P', p_4, "O2")
#     # Efficiency
#     eta_en = ( -(h_4-h_3)*(eta_mec_t*eta_pi) - (h_2-h_1)/(eta_mec_c*eta_pi) ) / ((h_3-h_2))
#
#     print('     | p\t| T \t\t| s-s1\t\t| h-h1')
#     print('     | [kPa]\t| [K] \t\t| [kJ/kg.K]\t| [kJ/kg]')
#     print(75*'-')
#     print('  1  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_1/1000, T_1, s_1/1000, h_1/1000))
#     print('  2  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_2/1000, T_2, s_2/1000, h_2/1000))
#     print('  3  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_3/1000, T_3, s_3/1000, h_3/1000))
#     print('  4  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_4/1000, T_4, s_4/1000, h_4/1000))
#     print(75*'-')
#     print('eta_en : %.2f' %eta_en)
#     print('\n')
#
#     # Final outputs - do not modify
#     p = (p_1, p_2, p_3, p_4)
#     T = (T_1, T_2, T_3, T_4)
#     s = (s_1, s_2, s_3, s_4)
#     h = (h_1, h_2, h_3, h_4)
#     out = (p,T,s,h,eta_en)
#     return out

def gas_turbine(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t):


    p_1, p_2, p_3, p_4 = p_1, p_2, p_3, p_1  # [Pa]
    T_1, T_2, T_3, T_4 = T_1, 0, T_3, 0    # [K]
    s_1, s_2, s_3, s_4 = 0, 0, 0, 0        # [J/kgK]
    h_1, h_2, h_3, h_4 = 0, 0, 0, 0        # [J/kg]

    eta_en = 0

    cp = 1.005e3          # [J/kgK]
    R_star = 8.3145e-3*(.21/32 + .79/28)
    # State 2
    T_2 = T_1*(p_2/p_1)**(R_star*cp/eta_pi)
    h_2 = cp*(T_2-T_1)
    # State 3
    h_3 = h_2 + cp*(T_3-T_2)
    s_3 = s_2 + cp*np.log(T_3/T_2)
    # State 4
    s_4 = s_3
    T_4 = T_3*(p_4/p_3)**(eta_pi*R_star/cp)
    h_4 = h_3 + cp*(T_4-T_3)

    # Efficiency
    eta_en = ( -(h_4-h_3)*eta_mec_t - (h_2-h_1)/eta_mec_c )*eta_pi / ((h_3-h_2))

    print('     | p\t| T \t\t| s-s1\t\t| h-h1')
    print('     | [kPa]\t| [K] \t\t| [kJ/kg.K]\t| [kJ/kg]')
    print(75*'-')
    print('  1  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_1/1000, T_1, s_1/1000, h_1/1000))
    print('  2  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_2/1000, T_2, s_2/1000, h_2/1000))
    print('  3  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_3/1000, T_3, s_3/1000, h_3/1000))
    print('  4  | %.2f\t| %.2f \t| %.2f \t\t| %.2f ' %(p_4/1000, T_4, s_4/1000, h_4/1000))
    print(75*'-')
    print('eta_en : %.2f' %eta_en)
    print('\n')

    # Final outputs - do not modify
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    s = (s_1, s_2, s_3, s_4)
    h = (h_1, h_2, h_3, h_4)
    out = (p,T,s,h,eta_en)
    return out

##########################################
p_1, T_1 = 1e+5, 293.15 # [Pa], [K]
p_2 = 17.8e+5 # [Pa]
p_3, T_3 = 17.8e+5, 1273.15 # [Pa], [K]
eta_pi = 0.90 # [-]
eta_mec_c, eta_mec_t = 0.98, 0.98 # [-], [-]
gas_turbine(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t)
##########################################

#
#===RANKINE CYCLE - TO BE IMPLEMENTED==========================================
#

def steam_turbine(T_1,p_3,T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump):
    # Replace with your model
    p_1, p_2, p_3, p_4 = 0, 0, 0, 0
    T_1, T_2, T_3, T_4 = 0, 0, 0, 0
    s_1, s_2, s_3, s_4 = 0, 0, 0, 0
    h_1, h_2, h_3, h_4 = 0, 0, 0, 0
    x_1, x_2, x_3, x_4 = 0, 0, 0, 0
    eta_en = 0
    dot_m_f = 0
    # Final outputs - do not modify
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    s = (s_1, s_2, s_3, s_4)
    h = (h_1, h_2, h_3, h_4)
    x = (x_1, x_2, x_3, x_4)
    out = (p,T,s,h,x,eta_en,dot_m_f)
    return out
