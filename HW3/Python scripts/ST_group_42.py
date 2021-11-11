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
#===CHECK BEFORE SUBMISSION====================================================
#   - clean the x_i == -1
#   - why is p_5 == 6200 [kPa] in the book (p91) and not 7000? -> 6200 may be hardcoded! -> don't forget to implement it correctly
#   - ditto for p_2 and p_3
#   - setting up the reference state; book seems to use the same one as the CoolProp library
#   - delete the initialization section
#   - check implementation of min/max value given in options
#   - coefficient 1.43 for p_1 arbitrarily chosen to get a pressure similar to what is obtained in the book

def print_red(str):
    print('\033[31m %s \033[0m' %str)
def print_green(str):
    print('\033[32m %s \033[0m' %str)

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
    p_1,p_2,p_3,p_4,p_5,p_6,p_6I,p_6II,p_6III,p_6IV,p_6V,p_6VI,p_6VII,p_6VIII, p_7,p_7I,p_7II,p_7III,p_7IV,p_7V,p_7VI,p_7VII,p_7VIII,p_8,p_9,p_9I,p_9II, p_9III,p_9IV,p_9V,p_9VI,p_9VII,p_9VIII = np.ones(9+8*3)*np.nan

    T_1,T_2,T_3,T_4,T_5,T_6,T_6I,T_6II,T_6III,T_6IV,T_6V,T_6VI,T_6VII,T_6VIII, T_7,T_7I,T_7II,T_7III,T_7IV,T_7V,T_7VI,T_7VII,T_7VIII,T_8,T_9,T_9I,T_9II, T_9III,T_9IV,T_9V,T_9VI,T_9VII,T_9VIII = np.ones(9+8*3)*np.nan

    s_1,s_2,s_3,s_4,s_5,s_6,s_6I,s_6II,s_6III,s_6IV,s_6V,s_6VI,s_6VII,s_6VIII, s_7,s_7I,s_7II,s_7III,s_7IV,s_7V,s_7VI,s_7VII,s_7VIII,s_8,s_9,s_9I,s_9II, s_9III,s_9IV,s_9V,s_9VI,s_9VII,s_9VIII = np.ones(9+8*3)*np.nan

    h_1,h_2,h_3,h_4,h_5,h_6,h_6I,h_6II,h_6III,h_6IV,h_6V,h_6VI,h_6VII,h_6VIII, h_7,h_7I,h_7II,h_7III,h_7IV,h_7V,h_7VI,h_7VII,h_7VIII,h_8,h_9,h_9I,h_9II, h_9III,h_9IV,h_9V,h_9VI,h_9VII,h_9VIII = np.ones(9+8*3)*np.nan

    e_1,e_2,e_3,e_4,e_5,e_6,e_6I,e_6II,e_6III,e_6IV,e_6V,e_6VI,e_6VII,e_6VIII, e_7,e_7I,e_7II,e_7III,e_7IV,e_7V,e_7VI,e_7VII,e_7VIII,e_8,e_9,e_9I,e_9II, e_9III,e_9IV,e_9V,e_9VI,e_9VII,e_9VIII = np.ones(9+8*3)*np.nan

    x_1,x_2,x_3,x_4,x_5,x_6,x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII, x_7,x_7I,x_7II,x_7III,x_7IV,x_7V,x_7VI,x_7VII,x_7VIII,x_8,x_9,x_9I,x_9II, x_9III,x_9IV,x_9V,x_9VI,x_9VII,x_9VIII = np.ones(9+8*3)*np.nan

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

    #===Ref. State=============================================================

    T_ref = T_ref   #[K]
    p_ref = p_ref   #[Pa]
    Dmolar = CP.PropsSI("Dmolar", "T", T_ref, "P", p_ref, 'Water')
    CP.set_reference_state('Water', T_ref, Dmolar, 0, 0)

    h_ref = CP.PropsSI('H','P',p_ref,'T', T_ref,'Water') #[J/kg]
    s_ref = CP.PropsSI('S','P',p_ref,'T',T_ref,'Water') #[J/kgK]
    e_ref = 0  #[J/kg]
    # print(' %.f\t%.3f\t%.3f' %(h_ref, s_ref, e_ref))
    #

    def exergy(h,s):
        return( (h-h_ref) - T_ref*(s-s_ref) )

    #
    #===CYCLE STATES===========================================================
    #

    # State 3 -- BOILER OUTPUT -- Transformations in boiler are supposed isobaric
    p_3 = p_3
    T_3 = T_max
    h_3 = CP.PropsSI('H','P',p_3,'T',T_3,'Water')
    s_3 = CP.PropsSI('S','P',p_3,'T',T_3,'Water')
    x_3 = CP.PropsSI('Q','P',p_3,'T',T_3,'Water')
    e_3 = exergy(h_3,s_3)

    # State 4 & 6VIII -- HIGH PRESSURE TURBINE OUTPUT -- Transformations in turbines are supposed adiabatic
    p_4 = p_4
    s_4is = s_3
    h_4is = CP.PropsSI('H','P',p_4,'S',s_4is,'Water')
    h_4 = h_3 - eta_is_HP*(h_3-h_4is)
    s_4 = CP.PropsSI('S','P',p_4,'H',h_4,'Water')
    T_4 = CP.PropsSI('T','P',p_4,'H',h_4,'Water')
    x_4 = CP.PropsSI('Q','P',p_4,'H',h_4,'Water')
    e_4 = exergy(h_4,s_4)

    # State 5 -- REHEATING -- cfr state 3
    # p_5 = 6200000
    p_5 = p_4
    T_5 = T_max
    h_5 = CP.PropsSI('H','P',p_5,'T',T_5,'Water')
    s_5 = CP.PropsSI('S','P',p_5,'T',T_5,'Water')
    x_5 = CP.PropsSI('Q','P',p_5,'T',T_5,'Water')
    e_5 = exergy(h_5,s_5)

    # State 6 -- CONDENSER INPUT
    x_6min = x_6
    T_6 = T_cond_out + T_pinch_cond
    p_6 = CP.PropsSI('P','T',T_6,'Q',1,'Water')
    s_6is = s_5
    h_6is = CP.PropsSI('H','P',p_6,'S',s_6is,'Water')
    h_6 = h_5 + (h_6is-h_5)*eta_is_LP
    x_6 = CP.PropsSI('Q','P',p_6,'H',h_6,'Water')
    s_6 = CP.PropsSI('S','P',p_6,'H',h_6,'Water')
    e_6 = exergy(h_6,s_6)

    if x_6 < x_6min:
        print_red('Warning: The steam fraction at the end of the turbine is lower than the minimum acceptable value:\n x_6 = %.2f\n Minimum acceptable value: %.2f' %(x_6, x_6min))
    else:
        print_green('The steam fraction at the end of the turbine is higher than the minimum acceptable value:\n x_6 = %.2f\n Minimum acceptable value: %.2f' %(x_6, x_6min))

    ## State 6 -- BLEEDS
    h_6_bleeds = np.linspace(h_5,h_6,9)[1:-1][::-1]
    p_6_bleeds = np.ones(7)*np.nan
    T_6_bleeds = np.ones(7)*np.nan
    x_6_bleeds = np.ones(7)*np.nan
    s_6_bleeds = np.ones(7)*np.nan
    e_6_bleeds = np.ones(7)*np.nan

    for i in np.arange(len(h_6_bleeds)):
        s_6is = s_5
        h_6is = h_5 + (h_6_bleeds[i]-h_5)/eta_is_LP
        p_6_bleeds[i] = CP.PropsSI('P','S',s_6is,'H',h_6is,'Water')
        T_6_bleeds[i] = CP.PropsSI('T','P',p_6_bleeds[i],'H',h_6_bleeds[i],'Water')
        x_6_bleeds[i] = CP.PropsSI('Q','P',p_6_bleeds[i],'H',h_6_bleeds[i],'Water')
        s_6_bleeds[i] = CP.PropsSI('S','P',p_6_bleeds[i],'H',h_6_bleeds[i],'Water')
        e_6_bleeds[i] = exergy(h_6_bleeds[i],s_6_bleeds[i])

    p_6_bleeds = np.append(p_6_bleeds,p_4)
    T_6_bleeds = np.append(T_6_bleeds,T_4)
    x_6_bleeds = np.append(x_6_bleeds,x_4)
    h_6_bleeds = np.append(h_6_bleeds,h_4)
    s_6_bleeds = np.append(s_6_bleeds,s_4)
    e_6_bleeds = np.append(e_6_bleeds,e_4)

    p_6I,p_6II,p_6III,p_6IV,p_6V,p_6VI,p_6VII,p_6VIII = p_6_bleeds
    T_6I,T_6II,T_6III,T_6IV,T_6V,T_6VI,T_6VII,T_6VIII = T_6_bleeds
    x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII = x_6_bleeds
    h_6I,h_6II,h_6III,h_6IV,h_6V,h_6VI,h_6VII,h_6VIII = h_6_bleeds
    s_6I,s_6II,s_6III,s_6IV,s_6V,s_6VI,s_6VII,s_6VIII = s_6_bleeds
    e_6I,e_6II,e_6III,e_6IV,e_6V,e_6VI,e_6VII,e_6VIII = e_6_bleeds

    # State 7 -- CONDENSER OUTPUT
    T_7 = T_cond_out
    x_7 = 0
    p_7 = CP.PropsSI('P','T',T_7,'Q',x_7,'Water')
    h_7 = CP.PropsSI('H','T',T_7,'Q',x_7,'Water')
    s_7 = CP.PropsSI('S','T',T_7,'Q',x_7,'Water')
    e_7 = exergy(h_7,s_7)

    ## State 7 -- HEAT EXCH.
    p_7_xch = p_6_bleeds
    T_7_xch = np.ones(8)*np.nan
    x_7_xch = np.zeros(8)
    h_7_xch = np.ones(8)*np.nan
    s_7_xch = np.ones(8)*np.nan
    e_7_xch = np.ones(8)*np.nan

    for i in np.arange(len(h_7_xch)):
        T_7_xch[i] = CP.PropsSI('T','P',p_7_xch[i],'Q',x_7_xch[i],'Water')
        h_7_xch[i] = CP.PropsSI('H','P',p_7_xch[i],'Q',x_7_xch[i],'Water')
        s_7_xch[i] = CP.PropsSI('S','P',p_7_xch[i],'Q',x_7_xch[i],'Water')
        e_7_xch[i] = exergy(h_7_xch[i],s_7_xch[i])

    p_7I,p_7II,p_7III,p_7IV,p_7V,p_7VI,p_7VII,p_7VIII = p_7_xch
    T_7I,T_7II,T_7III,T_7IV,T_7V,T_7VI,T_7VII,T_7VIII = T_7_xch
    x_7I,x_7II,x_7III,x_7IV,x_7V,x_7VI,x_7VII,x_7VIII = x_7_xch
    h_7I,h_7II,h_7III,h_7IV,h_7V,h_7VI,h_7VII,h_7VIII = h_7_xch
    s_7I,s_7II,s_7III,s_7IV,s_7V,s_7VI,s_7VII,s_7VIII = s_7_xch
    e_7I,e_7II,e_7III,e_7IV,e_7V,e_7VI,e_7VII,e_7VIII = e_7_xch

    if T_7IV < T_drum:
        print_red('Warning: The temperature inside de degassing drum is lower than the minimum acceptable value:\n T_7IV = %.2f\n Minimum acceptable value: %.2f' %(T_7IV, T_drum))
    else:
        print_green('The temperature inside de degassing drum is higher than the minimum acceptable value:\n T_7IV = %.2f [°C]\n Minimum acceptable value: %.2f [°C]' %(T_7IV-273.15, T_drum-273.15))

    # State 8 -- SECONDARY PUMP (P_e) OUTPUT -- Transformations in pumps are supposed adiabatic
    p_8 = p_7IV
    s_8is = s_7
    h_8is = CP.PropsSI('H','P',p_8,'S',s_8is,'Water')
    h_8 = h_7 + (h_8is-h_7)/eta_pump
    s_8 = CP.PropsSI('S','P',p_8,'H',h_8,'Water')
    T_8 = CP.PropsSI('T','P',p_8,'H',h_8,'Water')
    x_8 = CP.PropsSI('Q','P',p_8,'H',h_8,'Water')
    e_8 = exergy(h_8,s_8)

    # State 1 -- BOILER INPUT -- Transformations in heat ex. are supposed isobaric
    T_1 = T_7VIII-T_pinch_ex
    p_1 = CP.PropsSI('P','T',T_1,'Q',0,'Water')*1.43
    h_1 = CP.PropsSI('H','P',p_1,'T',T_1,'Water')
    s_1 = CP.PropsSI('S','P',p_1,'T',T_1,'Water')
    x_1 = CP.PropsSI('Q','P',p_1,'T',T_1,'Water')
    e_1 = exergy(h_1,s_1)

    # State 90 -- SUBCOOLER OUTPUT -- Transformations in heat ex. are supposed isobaric
    p_9 = p_8
    T_9 = T_8 + T_pinch_sub
    h_9 = CP.PropsSI('H','P',p_9,'T',T_9,'Water')
    s_9 = CP.PropsSI('S','P',p_9,'T',T_9,'Water')
    x_9 = CP.PropsSI('Q','P',p_9,'T',T_9,'Water')
    e_9 = exergy(h_9,s_9)

    ## State 9 -- HEAT EXCH.
    p_9_xch = np.ones(8)*np.nan
    T_9_xch = np.ones(8)*np.nan
    x_9_xch = np.zeros(8)
    h_9_xch = np.ones(8)*np.nan
    s_9_xch = np.ones(8)*np.nan
    e_9_xch = np.ones(8)*np.nan

    p_9_xch[0:3]  = np.ones( len(p_9_xch[0:3]) )*p_7IV
    p_9_xch[3:-1] = np.ones( len(p_9_xch[3:-1]))*p_1
    p_9_xch[-1] = p_1
    for i in np.arange(len(p_9_xch)):
        if i == 3:
            s_9IVis = s_7_xch[i]
            h_9IVis = CP.PropsSI('H','P',p_9_xch[i],'S',s_9IVis,'Water')
            h_9_xch[i] = h_7_xch[i] + (h_9IVis-h_7_xch[i])/eta_pump
            s_9_xch[i] = CP.PropsSI('S','P',p_9_xch[i],'H',h_9_xch[i],'Water')
            T_9_xch[i] = CP.PropsSI('T','P',p_9_xch[i],'H',h_9_xch[i],'Water')
            x_9_xch[i] = CP.PropsSI('Q','P',p_9_xch[i],'H',h_9_xch[i],'Water')
            e_9_xch[i] = exergy(h_9_xch[i],s_9_xch[i])
        else:
            T_9_xch[i] = T_7_xch[i]-T_pinch_ex
            x_9_xch[i] = CP.PropsSI('Q','P',p_9_xch[i],'T',T_9_xch[i],'Water')
            h_9_xch[i] = CP.PropsSI('H','P',p_9_xch[i],'T',T_9_xch[i],'Water')
            s_9_xch[i] = CP.PropsSI('S','P',p_9_xch[i],'T',T_9_xch[i],'Water')
            e_9_xch[i] = exergy(h_9_xch[i],s_9_xch[i])


    p_9I,p_9II,p_9III,p_9IV,p_9V,p_9VI,p_9VII,p_9VIII = p_9_xch
    T_9I,T_9II,T_9III,T_9IV,T_9V,T_9VI,T_9VII,T_9VIII = T_9_xch
    x_9I,x_9II,x_9III,x_9IV,x_9V,x_9VI,x_9VII,x_9VIII = x_9_xch
    h_9I,h_9II,h_9III,h_9IV,h_9V,h_9VI,h_9VII,h_9VIII = h_9_xch
    s_9I,s_9II,s_9III,s_9IV,s_9V,s_9VI,s_9VII,s_9VIII = s_9_xch
    e_9I,e_9II,e_9III,e_9IV,e_9V,e_9VI,e_9VII,e_9VIII = e_9_xch

    # State 2 -- BOILER INPUT -- Transformations in heat ex. are supposed isobaric
    #         -- MAIN PUMP OUTPUT -- Transformations in pumps are supposed adiabatic
    p_2 = p_3
    s_2is = s_1
    h_2is = CP.PropsSI('H','P',p_2,'S',s_2is,'Water')
    h_2 = h_1 + (h_2is-h_1)/eta_pump
    s_2 = CP.PropsSI('S','P',p_2,'H',h_2,'Water')
    T_2 = CP.PropsSI('T','P',p_2,'H',h_2,'Water')
    x_2 = CP.PropsSI('Q','P',p_2,'H',h_2,'Water')
    e_2 = exergy(h_2,s_2)
    
    
    ## Flow definition (Bleedings) 
    # ============================
    nsout=7
    reheat=1
    id_drum=3
    A = np.zeros((nsout+reheat,nsout+reheat));
    B = np.zeros(nsout+reheat);
    
    h9_b = np.array([h_9,h_9I,h_9II,h_9III,h_9IV,h_9V,h_9VI,h_9VII,h_9VIII])
    h7_b = np.array([h_7,h_7I,h_7II,h_7III,h_7IV,h_7V,h_7VI,h_7VII,h_7VIII])
    h6_b = np.array([h_6,h_6I,h_6II,h_6III,h_6IV,h_6V,h_6VI,h_6VII,h_6VIII])
    
    
    for i in range(id_drum):
        for j in range(id_drum):
            if (i==j):
                A[i][j] = (h9_b[i+1]-h9_b[i])-(h6_b[i+1]-h7_b[i+1]);
            elif (j>i):
                A[i][j] = (h9_b[i+1]-h9_b[i])-(h7_b[i+2]-h7_b[i+1]);
        B[i] = -(h9_b[i+1]-h9_b[i]);
    
    for j in range(nsout+reheat):
        if j < id_drum :
            A[id_drum][j] = (h7_b[id_drum+1] - h9_b[id_drum]);
        elif j==id_drum :
            A[id_drum][j] = -(h6_b[id_drum+1] - h7_b[id_drum]);
        else:
            A[id_drum][j] = -(h7_b[id_drum+1] - h7_b[id_drum+2]);
        B[id_drum] = -(h7_b[id_drum+1] - h9_b[id_drum]);
    for i in range(id_drum+1,nsout+reheat):
        for j in range(nsout+reheat):
            if (i==j):
                A[i][j] = (h9_b[i+1]-h9_b[i])-(h6_b[i+1]-h7_b[i+1]);
            elif (j>i):
                A[i][j] = (h9_b[i+1]-h9_b[i])-(h7_b[i+2]-h7_b[i+1]);
        B[i] = -(h9_b[i+1]-h9_b[i]);
    
    
    X = np.linalg.solve(A,B);
    Xtot = np.sum(X)
    
    
    #Turbine Work
    WmT = (1+Xtot)*(h_3-h_4)+(h_5-h_6);
    emT = (1+Xtot)*(e_3-e_4)+(e_5-e_6);
    
    #Pump Work
    WmP = (1+Xtot)*(h_2-h_1);
    emP = (1+Xtot)*(e_2-e_1);
    
    #Steam Generator
    Qh = (1+Xtot)*(h_3-h_2);
    eSG = (1+Xtot)*(e_3-e_2);
    
    #Condenser
    Qc = h_6-h_7;
    eCond = e_6-e_7;
    
    Qh += (1+Xtot-X[-1])*(h_5-h_4);
    eSG += (1+Xtot-X[-1])*(e_5-e_4);
    
    e9_b = np.array([e_9,e_9I,e_9II,e_9III,e_9IV,e_9V,e_9VI,e_9VII,e_9VIII])
    e7_b = np.array([e_7,e_7I,e_7II,e_7III,e_7IV,e_7V,e_7VI,e_7VII,e_7VIII])
    e6_b = np.array([e_6,e_6I,e_6II,e_6III,e_6IV,e_6V,e_6VI,e_6VII,e_6VIII])
            
    for i in range(1,nsout+reheat+1):
        WmT += X[i-1]*(h_5-h6_b[i]);
        emT += X[i-1]*(e_5-e6_b[i]);

        WmPe1 = (1+np.sum(X[:id_drum-1]))*(h_8-h_7);
        emPe1 = (1+np.sum(X[:id_drum-1]))*(e_8-e_7);
        WmPe2 = (1+Xtot)*(h9_b[id_drum]-h7_b[id_drum]);
        emPe2 = (1+Xtot)*(e9_b[id_drum]-e7_b[id_drum]);
        Qc += (1+np.sum(X[:id_drum-1]))*(h7_b[1]-h_7);
        eCond += (1+np.sum(X[:id_drum-1]))*(e7_b[1]-e_7);


    Wmtot = WmT-WmP-WmPe1-WmPe2 # Total work 
    eta_cyclen = Wmtot/Qh
    eta_cyclex = Wmtot/eSG

    # X_i COMPUTATION
    n = len(h_6_bleeds)
    A = np.zeros((n,n))
    b = np.zeros(n-1)

    for i in np.arange(n-1):
        b[i] = h_9_xch[i] - h_9_xch[i+1]
    b = np.append(np.append(b[:3], h_9_xch[3]-h_7_xch[3+1]),b[3:])
    # A[0][0] = (h_9_xch[1]-h_9_xch[0]) + (h_7_xch[0] - h_6_bleeds[1])
    # b[0] = h_9_xch[0]-h_9_xch[1]
    # for i in range(1,n-1):
    #     A[0][i]=(h9[1]-h9[0])+(h7[0]-h7[2])
    #     b[i] = h_9_xch[i]-h_9_xch[i+1]
    #     for j in range(1,n-1):
    #         if i > j:
    #             A[i][j] = (h_9_xch[i+1]-h_9_xch[i])
    #         elif i < j:
    #             A[i][j] = (h_9_xch[i+1]-h_9_xch[i])+(h_7_xch[i+1]-h_7_xch[i+2])
    #         elif i == j:
    #             A[i][j] = (h_9_xch[i+1]-h_9_xch[i])+(h_7_xch[i+1]-h_6_bleeds[i+1])

    print(A)
    print(75*'_')
    print(b, len(b))
    print(150*'_')


    #
    #===COMBUSTION=============================================================
    #




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
    x = (x_1,x_2,x_3,x_4,x_5,
          x_6,x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII,
          x_7,x_7I,x_7II,x_7III,x_7IV,x_7V,x_7VI,x_7VII,x_7VIII,
          x_8,
          x_9,x_9I,x_9II,x_9III,x_9IV,x_9V,x_9VI,x_9VII,x_9VIII)

    x = np.array(x)
    x[x==-1] = np.ones(len(x[x==-1]))*np.nan
    x = tuple(x)

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
