"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Signature of the final steam turbine function

@author: Antoine Laterre
@date: October 21, 2021
"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

#
#===STEAM CYCLE - TO BE IMPLEMENTED============================================
#

#
#===STEAM CYCLE - TO BE IMPLEMENTED============================================
#

# SPECIFICATIONS:
#   steam_turbine(P_e,options,display) computes the thermodynamics states for
#   the steam power plant of pp. 89, fig 2.33 (combustion, exchanger, cycle)
#   based on several inputs (given in options) and based on a given electricity
#   production P_e.  It returns the main results. It can also plot the T-s and
#   h-s as well as the energy and exergy pies if display == True
# INPUTS (/!\ some inputs can be dependent on others /!\):
#   - P_e               [W]: net power output
#   - options: tuple containing the parametric values of the GT
#       o p3            [Pa]: maximum steam pressure
#       o p4            [Pa]: reheating pressure
#       o p_ref         [Pa]: reference pressure for exergy computation
#       o T_ref         [K]: reference temperature for exergy computation
#       o T_max         [K]: maximum steam temperature
#       o T_cond_out    [K]: condenser cold outlet temperature
#       o T_exhaust     [K]: exhaust gas temperature out of the chimney
#       o T_pinch_sub   [K]: pinch temperature at the subcooler
#       o T_pinch_ex    [K]: pinch temperature at heat exchangers
#       o T_pinch_cond  [K]: pinch temperature at the condenser
#       o T_drum        [K]: minimum drum temperature
#       o x_6           [-]: minimum possible vapor quality after final expansion
#       o comb          [na]: tuple containing combustion data
#           * Tmax      [K]: maximum combustion temperature
#           * lambda    [-]: excess air coefficient
#           * x         [-]: O_x/C ratio of fuel (e.g. 0.05 in CH_1.2O_0.05)
#           * y         [-]: H_y/C ratio of fuel (e.g. 1.20 in CH_1.2O_0.05)
#       o eta_mec       [-]: shafts bearings mechanical efficiency
#       o eta_pump      [-]: internal efficiency of the pump (see text book pp. 53)
#       o eta_is_turb   [na]: tuple containing turbines isentropic efficiencies (see text book pp. 55)
#           * eta_is_HP [-]: high pressure turbine isentropic efficiency
#           * eta_is_LP [-]: low pressure turbine isentropic efficiency
#   - display: bool to choose to plot the T-s & h-s diagrams and the energy and exergy pie charts (True or False)
#
# OUTPUTS: tuple containing...
#   - DAT: tuple containing the cycle state data (1,2,3, ... , 6,6_I, ... ,9_VII,9_VIII)
#      o p              [Pa]: tuple containing the pressure at each state
#      o T              [K]: tuple containing the temperature at each state
#      o h              [J/kg]: tuple containing the enthalpy at each state
#      o s              [J/kg/K]: tuple containing the entropy at each state
#      o e              [J/kg]: tuple containing the exergy at each state (use ref conditions)
#      o x              [-]: tuple containing the vapor quality at each state
#   - COMBUSTION: tuple containing the combustion parameters
#      o LHV            [J/kg]: the fuel Lower Heating Value
#      o e_c            [J/kg]: the fuel exergy
#      o excess_air     [-]: excess air of the combustion
#      o cp_gas         [J/kg/K]: the flue gas specific heat capacity at the combustor outlet
#      o gas_prop       [-]: list containing the proportion of ['N2','O2','CO2','H2O'] in the flue gas respectively
#   - MASSFLOW: tuple containing the massflow rates
#      o dot_m_a        [kg/s]: mass flow rate of air
#      o dot_m_f        [kg/s]: mass flow rate of fuel
#      o dot_m_g        [kg/s]: mass flow rate of flue gas
#      o dot_m_v        [kg/s]: mass flow rate of steam
#   - ETA: tuple containing the efficiencies (see text book pp. 53-94)
#       o eta_cyclen    [-]: cycle energy efficiency
#       o eta_toten     [-]: overall energy efficiency
#       o eta_cyclex    [-]: cycle exergy efficiency
#       o eta_totex     [-]: overall exergy efficiency
#       o eta_gen       [-]: steam generator energy efficiency
#       o eta_gex       [-]: steam generator exergy efficiency
#       o eta_combex    [-]: combustion exergy efficiency
#       o eta_chimnex   [-]: chimney exergy efficiency
#       o eta_condex    [-]: condenser exergy efficiency
#       o eta_transex   [-]: bleedings heat exchangers overall exergy efficiency
#       o eta_rotex     [-]: pumps and turbines exergy efficiency
#   - DATEN: tuple containing the energy losses
#       o loss_mec      [W]: mechanical energy losses
#       o loss_gen      [W]: steam generator energy losses
#       o loss_cond     [W]: condenser energy losses
#   - DATEX: tuple containing the exergy losses
#       o loss_mec      [W]: mechanical energy losses
#       o loss_rotex    [W]: pumps and turbines exergy losses
#       o loss_combex   [W]: combustion exergy losses
#       o loss_chemex   [W]: chimney exergy losses
#       o loss_transex  [W]: bleedings heat exchangers overall exergy losses
#       o loss_totex    [W]: total exergy losses
#       o loss_condex   [W]: condenser exergy losses
#   - XMASSFLOW: tuple containing the massflow rates in each feedheater (w.r.t. fig. 2.33, pp. 91)
#       o x_6I          [kg/s]:
#       o x_6II         [kg/s]:
#       o x_6...        [kg/s]:
#       o x_6VIII       [kg/s]:
#   - FIG: tuple containing the figures to be diplayed
#       o fig_pie_en: pie chart of energy losses
#       o fig_pie_ex: pie chart of exergy losses
#       o fig_Ts_diagram: T-s diagram of the GT cycle
#       o fig_hs_diagram: h-s diagram of the GT cycle
def steam_turbine(P_e,options,display):

    ###### initialisation -- delete next section when code is functioning as intended #####
    ##### SECTION #####
    p_1,p_2,p_3,p_4,p_5,p_6,p_6I,p_6II,p_6III,p_6IV,p_6V,p_6VI,p_6VII,p_6VIII, p_7,p_7I,p_7II,p_7III,p_7IV,p_7V,p_7VI,p_7VII,p_7VIII,p_8,p_9,p_9I,p_9II, p_9III,p_9IV,p_9V,p_9VI,p_9VII,p_9VIII = np.zeros(9+8*3)

    T_1,T_2,T_3,T_4,T_5,T_6,T_6I,T_6II,T_6III,T_6IV,T_6V,T_6VI,T_6VII,T_6VIII, T_7,T_7I,T_7II,T_7III,T_7IV,T_7V,T_7VI,T_7VII,T_7VIII,T_8,T_9,T_9I,T_9II, T_9III,T_9IV,T_9V,T_9VI,T_9VII,T_9VIII = np.zeros(9+8*3)

    s_1,s_2,s_3,s_4,s_5,s_6,s_6I,s_6II,s_6III,s_6IV,s_6V,s_6VI,s_6VII,s_6VIII, s_7,s_7I,s_7II,s_7III,s_7IV,s_7V,s_7VI,s_7VII,s_7VIII,s_8,s_9,s_9I,s_9II, s_9III,s_9IV,s_9V,s_9VI,s_9VII,s_9VIII = np.zeros(9+8*3)

    h_1,h_2,h_3,h_4,h_5,h_6,h_6I,h_6II,h_6III,h_6IV,h_6V,h_6VI,h_6VII,h_6VIII, h_7,h_7I,h_7II,h_7III,h_7IV,h_7V,h_7VI,h_7VII,h_7VIII,h_8,h_9,h_9I,h_9II, h_9III,h_9IV,h_9V,h_9VI,h_9VII,h_9VIII = np.zeros(9+8*3)

    e_1,e_2,e_3,e_4,e_5,e_6,e_6I,e_6II,e_6III,e_6IV,e_6V,e_6VI,e_6VII,e_6VIII, e_7,e_7I,e_7II,e_7III,e_7IV,e_7V,e_7VI,e_7VII,e_7VIII,e_8,e_9,e_9I,e_9II, e_9III,e_9IV,e_9V,e_9VI,e_9VII,e_9VIII = np.zeros(9+8*3)

    x_1,x_2,x_3,x_4,x_5,x_6,x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII, x_7,x_7I,x_7II,x_7III,x_7IV,x_7V,x_7VI,x_7VII,x_7VIII,x_8,x_9,x_9I,x_9II, x_9III,x_9IV,x_9V,x_9VI,x_9VII,x_9VIII = np.zeros(9+8*3)

    LHV,e_c,excess_air,cp_gas,gas_prop = np.zeros(5)
    dotm_a,dotm_f,dotm_g,dotm_v = np.zeros(4)
    eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex = np.zeros(6)
    loss_mec,loss_gen,loss_cond = np.zeros(3)
    loss_mec,loss_rotex,loss_combex,loss_chemex,loss_transex,loss_totex,loss_condex = np.zeros(7)
    fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram = [],[],[],[]
    ###### END OF SECTION #####

    # Process input variables--------------------------------------------------
    p_3,p_4,p_ref,T_ref,T_max,T_cond_out,T_exhaust,T_pinch_sub,T_pinch_ex,T_pinch_cond,T_drum,x_6,comb,eta_mec,eta_pump,eta_turb = options
    Tmax,exc_air,x,y = comb
    eta_is_HP,eta_is_LP = eta_turb

    # Ref. State
    T_ref = T_ref   #[K]
    p_ref = p_ref   #[Pa]
    Dmolar = CP.PropsSI("Dmolar", "T", T_ref, "P", p_ref, 'Water')
    CP.set_reference_state('Water', T_ref, Dmolar, 0, 0)

    h_ref = CP.PropsSI('H','P',p_ref,'T', T_ref,'Water') #[J/kg]
    s_ref = CP.PropsSI('S','P',p_ref,'T',T_ref,'Water') #[J/kgK]
    e_ref = 0; # Reference state #[J/kg]
    print(' %.f\t%.3f\t%.3f' %(h_ref, s_ref, e_ref))


    def exergy(h,s):
        return( (h-h_ref) - T_ref*(s-s_ref) )

    # CYCLE STATES

    # State 3 -- BOILER OUTPUT -- Transformations in boiler are supposed isobaric
    p_3 = p_3
    T_3 = T_max
    h_3 = CP.PropsSI('H','P',p_3,'T',T_3,'Water')
    s_3 = CP.PropsSI('S','P',p_3,'T',T_3,'Water')
    x_3 = CP.PropsSI('Q','P',p_3,'T',T_3,'Water')
    e_3 = exergy(h_3,s_3)

    # State 4 -- HIGH PRESSURE TURBINE OUTPUT -- Transformations in turbine are supposed adiabatic
    p_4 = p_4
    s_4is = s_3
    h_4is = CP.PropsSI('H','P',p_4,'S',s_4is,'Water')
    h_4 = h_3 - eta_is_HP*(h_3-h_4is)
    s_4 = CP.PropsSI('S','P',p_4,'H',h_4,'Water')
    T_4 = CP.PropsSI('T','P',p_4,'H',h_4,'Water')
    x_4 = CP.PropsSI('Q','P',p_4,'H',h_4,'Water')
    e_4 = exergy(h_4,s_4)

    # State 5 -- REHEATING -- cfr state 3
    p_5 = p_4
    T_5 = T_max                                         # why 6200 in the book (p91) and not 7000?
    h_5 = CP.PropsSI('H','P',p_5,'T',T_5,'Water')
    s_5 = CP.PropsSI('S','P',p_5,'T',T_5,'Water')
    x_5 = CP.PropsSI('Q','P',p_5,'T',T_5,'Water')
    e_5 = exergy(h_5,s_5)

    # State 6 -- CONDENSER INPUT
    T_6 = T_cond_out + T_pinch_cond
    x_6 = x_6
    p_6 = CP.PropsSI('P','T',T_6,'Q',x_6,'Water')
    h_6 = CP.PropsSI('H','T',T_6,'Q',x_6,'Water')
    s_6 = CP.PropsSI('S','T',T_6,'Q',x_6,'Water')
    e_6 = exergy(h_6,s_6)

    # State 7 -- CONDENSER OUTPUT
    T_7 = T_cond_out
    x_7 = 0
    p_7 = CP.PropsSI('P','T',T_7,'Q',x_7,'Water')
    h_7 = CP.PropsSI('H','T',T_7,'Q',x_7,'Water')
    s_7 = CP.PropsSI('S','T',T_7,'Q',x_7,'Water')
    e_7 = exergy(h_7,s_7)




    # State 1 -- MAIN PUMP INPUT




    # Process output variables - do not modify---------------------------------
    p = (p_1,p_2,p_3,p_4,p_5,
          p_6,p_6I,p_6II,p_6III,p_6IV,p_6V,p_6VI,p_6VII,p_6VIII,
          p_7,p_7I,p_7II,p_7III,p_7IV,p_7V,p_7VI,p_7VII,p_7VIII,
          p_8,
          p_9,p_9I,p_9II,p_9III,p_9IV,p_9V,p_9VI,p_9VII,p_9VIII)
    T = (T_1,T_2,T_3,T_4,T_5,
          T_6,T_6I,T_6II,T_6III,T_6IV,T_6V,T_6VI,T_6VII,T_6VIII,
          T_7,T_7I,T_7II,T_7III,T_7IV,T_7V,T_7VI,T_7VII,T_7VIII,
          T_8,
          T_9,T_9I,T_9II,T_9III,T_9IV,T_9V,T_9VI,T_9VII,T_9VIII)
    s = (s_1,s_2,s_3,s_4,s_5,
          s_6,s_6I,s_6II,s_6III,s_6IV,s_6V,s_6VI,s_6VII,s_6VIII,
          s_7,s_7I,s_7II,s_7III,s_7IV,s_7V,s_7VI,s_7VII,s_7VIII,
          s_8,
          s_9,s_9I,s_9II,s_9III,s_9IV,s_9V,s_9VI,s_9VII,s_9VIII)
    h = (h_1,h_2,h_3,h_4,h_5,
          h_6,h_6I,h_6II,h_6III,h_6IV,h_6V,h_6VI,h_6VII,h_6VIII,
          h_7,h_7I,h_7II,h_7III,h_7IV,h_7V,h_7VI,h_7VII,h_7VIII,
          h_8,
          h_9,h_9I,h_9II,h_9III,h_9IV,h_9V,h_9VI,h_9VII,h_9VIII)
    e = (e_1,e_2,e_3,e_4,e_5,
          e_6,e_6I,e_6II,e_6III,e_6IV,e_6V,e_6VI,e_6VII,e_6VIII,
          e_7,e_7I,e_7II,e_7III,e_7IV,e_7V,e_7VI,e_7VII,e_7VIII,
          e_8,
          e_9,e_9I,e_9II,e_9III,e_9IV,e_9V,e_9VI,e_9VII,e_9VIII)
    x = (x_1,x_2,x_3,x_4,x_5,                                           # PAS OUBLIER DE CONVERTIR LES x=-1 en NaN ou "-"
          x_6,x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII,
          x_7,x_7I,x_7II,x_7III,x_7IV,x_7V,x_7VI,x_7VII,x_7VIII,
          x_8,
          x_9,x_9I,x_9II,x_9III,x_9IV,x_9V,x_9VI,x_9VII,x_9VIII)
    DAT = (p,T,h,s,e,x)
    COMBUSTION = (LHV,e_c,excess_air,cp_gas,gas_prop)
    MASSFLOW = (dotm_a,dotm_f,dotm_g,dotm_v)
    ETA = (eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex)
    DATEN = (loss_mec,loss_gen,loss_cond)
    DATEX = (loss_mec,loss_rotex,loss_combex,loss_chemex,loss_transex,loss_totex,loss_condex)
    XMASSFLOW = (x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII)
    FIG = (fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram)
    out = (ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG)
    return out
