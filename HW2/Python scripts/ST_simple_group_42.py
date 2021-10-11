"""
LELME2150 - Thermal cycles
Homework 2 - Exergy and combustion

Signature of the simple steam turbine function

@author: Antoine Laterre
@date: September 22, 2021
"""

#
#===IMPORT PACKAGES============================================================
#

import CoolProp.CoolProp as CP
import numpy as np


#
#===RANKINE CYCLE - TO BE IMPLEMENTED==========================================
#

# IN:
#   - T_1 [K]: condenser hot leg outlet/pump inlet temperature
#   - p_3 [Pa], T_3 [K]: steam generator outlet pressure and temperature
#   - eta_gen [-]: steam generator energy efficiency (see text book pp. 57)
#   - eta_mec [-]: mechanical efficiency of the system (see text book pp. 58)
#   - eta_is_t [-]: isentropic efficiency of the turbine (see text book pp. 55)
#   - eta_pump [-]: internal efficiency of the pump (see text book pp. 53)
#   - LHV [J/kg]: Lower Heating Value of the fuel
#   - P_e [W]: net power output of the steam turbine
# OUT: tuple containing
#   - p [Pa]: tuple containing the pressure at each state
#   - T [T]: tuple containing the temperature at each state
#   - h [J/kg]: tuple containing the enthalpy at each state
#   - s [J/kg/K]: tuple containing the entropy at each state
#   - e [J/kg]: tuple containing the exergy at each state
#   - x [-]: tuple containing the vapor quality at each state ("titre" in French)
#   - eta_en [-]: tuple containing the energy efficiencies
#       o eta_cyclen [-]: cycle energy efficiency (see text book pp. 58)
#       o eta_toten [-]: overall energy efficiency (see text book pp. 5)
#   - dot_m_f [kg/s]: the fuel flow rate
# REMARK:
#   - You can adapt the numerical application of the text book in page pp. 62
#   to check your results.
def steam_turbine(T_1,p_3,T_3,eta_gen,eta_mec,eta_is_t,eta_pump,LHV,P_e):
    
    # Replace with your model--------------------------------------------------
    #p_1, p_2, p_3, p_4 = 0, 0, 0, 0
    #T_1, T_2, T_3, T_4 = 0, 0, 0, 0
    #h_1, h_2, h_3, h_4 = 0, 0, 0, 0
    #s_1, s_2, s_3, s_4 = 0, 0, 0, 0
    #e_1, e_2, e_3, e_4 = 0, 0, 0, 0
    #x_1, x_2, x_3, x_4 = 0, 0, 0, 0
    #eta_cyclen, eta_toten = 0, 0
    #eta_en = (eta_cyclen, eta_toten)
    #dot_m_f = 0
    
    #State 1: Saturated liquid
    T_1 = T_1
    x_1 = 0
    p_1 = CP.PropsSI('P','T',T_1,'Q',0,'Water')
    h_1 = CP.PropsSI('H','P',p_1,'Q',0,'Water')
    s_1 = CP.PropsSI('S','P',p_1,'Q',0,'Water')
    e_1 = 0; # Reference state

    #State 2: Sub-cooled water (1->2 supposed isothermal, 2->3 supposed isobaric)
    p_2 = p_3
    rho = CP.PropsSI('D','T',T_1,'P',p_1,'Water')
    T_2 = T_1 + (p_2-p_1)*(1/eta_pump -1)/(rho*CP.PropsSI('C','P',p_1,'T',T_1,'Water'))
    h_2 = CP.PropsSI('H','P',p_2,'T',T_2,'Water')
    s_2 = CP.PropsSI('S','P',p_2,'T',T_2,'Water')
    x_2 = CP.PropsSI('Q','P',p_2,'T',T_2,'Water')
    e_2 = (h_2 - h_1) - T_1*(s_2 - s_1)

    #State 3:
    p_3 = p_3
    T_3 = T_3
    h_3 = CP.PropsSI('H','P',p_3,'T',T_3,'Water')
    s_3 = CP.PropsSI('S','P',p_3,'T',T_3,'Water')
    x_3 = CP.PropsSI('Q','P',p_3,'T',T_3,'Water')
    e_3 = (h_3 - h_1) - T_1*(s_3 - s_1)

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
    e_4 = (h_4 - h_1) - T_1*(s_4 - s_1)
    
    #Mechanical work of the turbine:
    wm_t = h_3 - h_4
    
    #Mechanical work of the pump:
    wm_p = h_2-h_1
    
    #Heat
    q = h_3-h_2

    #Overall efficiency of the cycle:
    eta_cyclen = (wm_t-wm_p)/q

    eta_toten = eta_mec*eta_cyclen*eta_gen;
    
    eta_en = (eta_cyclen, eta_toten)

    #Water flow in the cycle:
    dot_m_v = P_e/(eta_mec*(wm_t-wm_p))

    #Mass flow of coal used:
    dot_m_f = dot_m_v*q/(eta_gen*LHV)
    
    #Losses:
    
    L_mec = dot_m_v*(wm_t-wm_p) - P_e; # Mechanical losses [W]
    
    L_cond = dot_m_v*(h_4-h_1); # Condenser losses [W]
    
    L_gen = dot_m_f*LHV - dot_m_v*q; # Steam generator losses [W]
    
    # Final outputs - do not modify--------------------------------------------
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    h = (h_1, h_2, h_3, h_4)
    s = (s_1, s_2, s_3, s_4)
    e = (e_1, e_2, e_3, e_4)
    x = (x_1, x_2, x_3, x_4)
    out = (p,T,h,s,e,x,eta_en,dot_m_f)
    return out
