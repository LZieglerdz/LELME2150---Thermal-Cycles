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

# Import yours (CoolProp, ...)

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
    p_1, p_2, p_3, p_4 = 0, 0, 0, 0
    T_1, T_2, T_3, T_4 = 0, 0, 0, 0
    h_1, h_2, h_3, h_4 = 0, 0, 0, 0
    s_1, s_2, s_3, s_4 = 0, 0, 0, 0
    e_1, e_2, e_3, e_4 = 0, 0, 0, 0
    x_1, x_2, x_3, x_4 = 0, 0, 0, 0
    eta_cyclen, eta_toten = 0, 0
    eta_en = (eta_cyclen, eta_toten)
    dot_m_f = 0
    
    # Final outputs - do not modify--------------------------------------------
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    h = (h_1, h_2, h_3, h_4)
    s = (s_1, s_2, s_3, s_4)
    e = (e_1, e_2, e_3, e_4)
    x = (x_1, x_2, x_3, x_4)
    out = (p,T,h,s,e,x,eta_en,dot_m_f)
    return out