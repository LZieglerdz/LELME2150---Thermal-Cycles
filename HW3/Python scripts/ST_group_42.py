"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Signature of the final steam turbine function

@author: Antoine Laterre
@date: October 21, 2021
"""

import CoolProp.CoolProp as CP
import numpy as np
from thermochem import janaf
db = janaf.Janafdb();
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
#   - check output units


#
#===GLOBAL VARIABLES===========================================================
#
T_S = 273.153         # [K]
p_S = 1e5             # [Pa]

# air composition
N2_conc = .79         # %mol
O2_conc = .21         # %mol
H2O_conc = 0          # %mol
CO2_conc = 0          # %mol

Mm_O2 = CP.PropsSI("MOLARMASS", "O2")        # kg/mol
Mm_N2 = CP.PropsSI("MOLARMASS", "N2")        # kg/mol
Mm_CO2 = CP.PropsSI("MOLARMASS", "CO2")        # kg/mol
Mm_H2O = CP.PropsSI("MOLARMASS", "H2O")        # kg/mol
Mm_H2 = CP.PropsSI("MOLARMASS", "H2")        # kg/mol
Mm_Ar = CP.PropsSI("MOLARMASS", "Ar")        # kg/mol
Mm_C = Mm_CO2-Mm_O2        # kg/mol
Mm_H = .5*Mm_H2         # kg/mol
Mm_O = .5*Mm_O2      # kg/mol

comp = ["N2","O2","CO2","H2O"]
air_conc = [N2_conc,O2_conc,CO2_conc,H2O_conc]
Mm   = [Mm_N2,Mm_O2,Mm_CO2,Mm_H2O]
Mm_air = np.dot(Mm, air_conc)

def print_red(str):
    print('\033[31m %s \033[0m' %str)
def print_green(str):
    print('\033[32m %s \033[0m' %str)



#
#===STEAM CYCLE================================================================
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
    p_5 = p_4
    p_5 = 6200000
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
    p_5 = 35000000
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


    #
    #===COMBUSTION=============================================================
    #

    def get_ma1(x, y): #Stoechiometric air-to-fuel ratio [-]
        return (Mm_O2+3.76*Mm_N2)*(1+(y-2*x)*0.25)/(Mm_C+Mm_H*y+Mm_O*x);

    def get_LHV(x, y): # Lower Heating Value for solid fuels CHyOx
        LHV_mol = (393400 + 102250*y - x*(111000 + 102250*y)/(1 + y/2))*1e3 # [J/kmol]
        LHV = LHV_mol/(12 + y + 16*x) # [J/kg]
        return LHV

    def CombEx(x, y): # CHyOx # e_c and Cp of different fuels (from slides of LMECA2150)
        if x==0: # Combustible of type CHy
            if y==0: # C
                e_c = 34160*1e3 # [J/kg]
                Cp = 0.86667*1e3
            if y==1.8: # CH1.8
                e_c = 45710*1e3 # [J/kg]
                Cp = 1.10434*1e3
            if y==4: # CH4
                e_c = 52215*1e3 # [J/kg]
                Cp = 2.2005*1e3
        elif y==0: # Combustible of type COx
            if x==1: # CO
                e_c = 9845*1e3 # [J/kg]
                Cp = 1.71176*1e3
        else:
            m_water = (y/2)*18/1000; # [kg]
            HHV = (get_LHV(x,y) + 2375*m_water*1e3) # [J/kg]
            e_c = HHV - 5350*1e3; # [J/kg] T0*S = 5350e3 (S is the carbon entropy at standard conditions)
            Cp = 0.86667*1e3
        return [Cp, e_c]

    def FlueGasFrac(x, y, lamb):
        w = lamb*(1+0.25*(y-2*x));
        x_O2f = (lamb-1)*(1+0.25*(y-2*x));
        x_N2f = 3.76*w;
        x_CO2f = 1;
        x_H2Of = y/2;
        x_f = [x_O2f,x_N2f,x_CO2f,x_H2Of];
        x_O2f = x_O2f/sum(x_f);
        x_N2f = x_N2f/sum(x_f);
        x_CO2f = x_CO2f/sum(x_f);
        x_H2Of = x_H2Of/sum(x_f);
        Mm_f = Mm_O2*x_O2f + Mm_N2*x_N2f + Mm_CO2*x_CO2f + Mm_H2O*x_H2Of; # Molar mass of fluegas [g/mol]
        R_f = janaf.R/Mm_f; # Ideal gas constant (R*) for exhaust gas [J/kg/K]

        return [Mm_f, [x_O2f, x_N2f, x_CO2f, x_H2Of], R_f]

    # def get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3,iter, lam_est):
    #     #  Combustion supposée complète
    #     #  C_zH_yO_x + w(O_2 + 3.76 N_2) --> a_0 O_2 + a_2 CO_2 + b_2 H_2O + 3.76w N_2
    #     # comb stoechiometric -> w = z + (y-2*x)/4
    #     # a_2 = z
    #     # 2*b_2 = y
    #     # 2*a_0 + 2*a_2 + b_2 = x + 2*w
    #     w = lam_est*(z + (y-2*x)/4)
    #     a_0 = (x + 2*w - (y/2) - 2*z)/2
    #     Mm_f, flue_conc_mol, R_f = FlueGasFrac(x, y, lam_est)
    #     flue_conc_mass = np.multiply(flue_conc_mol, [Mm_O2, Mm_N2, Mm_CO2, Mm_H2O])/Mm_f
    #
    #     cp_3S = getMeanCp(p_S, p_3, T_S, T_3, R_f, comp, flue_conc_mass)
    #     cp_32 = getMeanCp(p_2, p_3, T_2, T_3, R_f, comp, flue_conc_mass)
    #
    #     lam = (get_LHV(x,y)*1e3-cp_3S*(T_3-T_S))/ (get_ma1(x,y)*cp_32*(T_3-T_2))
    #
    #     if iter == 0:
    #         print('The function get_lambda did not converge.')
    #         return(-1)
    #     if np.abs(lam-lam_est) <= 1e-6 :
    #         return(lam)
    #     else:
    #         return( get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3, iter-1, lam) )

    def getCpMix(T, p, mix_comp, mix_conc):
        sum = 0
        for i in range(len(mix_comp)):
             sum += (mix_conc[i]*CP.PropsSI('CPMASS', 'T', T, 'P', p, mix_comp[i]) )
        return(sum)

    def getMeanCp(p_in, p_out, T_in, T_out, R, mix_comp, mix_conc):               #renvoie le cp massique moyen
        cp = getCpMix(T_in, p_in, mix_comp, mix_conc)

        if T_in == T_out:
            return(cp)
        n=100
        rangeT = np.linspace(T_in,T_out,n)
        rangeP = np.ones(n)*p_in
        p = p_in
        for i in np.arange(1,n):
            gamma = cp/(cp - i*R)
            p = p_in*(rangeT[i]/T_in)**((gamma-1)/gamma)
            cp += getCpMix(rangeT[i], rangeP[i], mix_comp, mix_conc)
        return(cp/n)

    ma1 = get_ma1(x,y)
    T_max_comb,lamb,x,y = comb
    Mm_f, [x_O2f, x_N2f, x_CO2f, x_H2Of], R_f = FlueGasFrac(x,y,lamb)

    lamb_ma1 = lamb*ma1

    LHV = get_LHV(x, y)                         #[J/kg]: the fuel Lower Heating Value
    # cp_gas                                    #[J/kg/K]: the flue gas specific heat capacity at the combustor outlet
    cp_gas, e_c  = CombEx(x, y)                 #[J/kg]: the fuel exergy
    excess_air =  lamb                          #[-]: excess air of the combustion
    gas_prop = [x_N2f, x_O2f, x_CO2f, x_H2Of]   #[-]: list containing the proportion of ['N2','O2','CO2','H2O'] in the flue gas respectively

    h_exh = getMeanCp(p_ref, p_ref, T_ref, T_exhaust, R_f, comp, gas_prop)*(T_exhaust-T_ref)#*1e-3 #[kJ/kg]
    h_a = getMeanCp(p_ref, p_ref, T_ref, T_ref, R_f, comp, air_conc)*(T_ref-T_ref)#*1e-3 #[kJ/kg]
    eps_exh = ((lamb_ma1+1)*h_exh - lamb_ma1*h_a)/LHV
    eps_p = .01                                 # "we assume there are no unburnt residues etc" -p63
    eta_gen = 1-eps_p-eps_exh
    print('eps_exh: %.2f [-]' %eps_exh)
    print('eps_p: %.2f [-]' %eps_p)
    print('eta_gen: %.2f [-]' %eta_gen)
    print('Qh: %.2f [J/kg]' %Qh)
    print('LHV: %.2f [J/kg]' %LHV)

    dot_m_vc = P_e / ( eta_mec * Wmtot )
    dot_m_vG = dot_m_vc*(1+Xtot)

    dot_m_v = dot_m_vc                          #[kg/s]: mass flow rate of steam
    dot_m_f = dot_m_vG*Qh/(eta_gen*LHV)         #[kg/s]: mass flow rate of fuel
    dot_m_a = dot_m_f*lamb_ma1                  #[kg/s]: mass flow rate of air
    dot_m_g = dot_m_a+dot_m_f                   #[kg/s]: mass flow rate of flue gas

    print('Steam massflow (cond.): %.2f [kg/s]' %dot_m_v)
    print('Steam massflow (boil.): %.2f [kg/s]' %dot_m_vG)
    print('Air massflow: %.2f [kg/s]' %dot_m_a)
    print('Fuel massflow: %.2f [kg/s]' %dot_m_f)
    print('Flue gas massflow: %.2f [kg/s]' %dot_m_g)
    print(75*'_')

    #
    #===EFFICIENCIES===========================================================
    #
    e_r = 0
    c_pf = getMeanCp(p_ref, p_ref, T_ref, T_max_comb, R_f, comp, gas_prop)
    c_pf_mean = c_pf/np.log(T_max_comb/T_ref)
    e_f = LHV/(lamb_ma1 + 1) - c_pf_mean*T_ref*np.log(1 + LHV/((lamb_ma1 + 1)*c_pf*T_ref ) )
    e_exh = 0
    for i in range(len(gas_prop)):
        e_exh += gas_prop[i]*( (CP.PropsSI('H','P',p_ref,'T',T_exhaust,comp[i]) - CP.PropsSI('H','P',p_ref,'T',T_ref,comp[i])) + T_ref*(CP.PropsSI('S','P',p_ref,'T',T_exhaust,comp[i])-CP.PropsSI('S','P',p_ref,'T',T_ref,comp[i]) ) )

    # print(25*'_')
    # print('e_f (flue gasses after combustion): %.3f [kJ/kg]' %(e_f*1e-3))
    # print('e_exh (flue gasses at chimney): %.3f [kJ/kg]' %(e_exh*1e-3))
    # print(25*'_')


    Cp_w = CP.PropsSI('C','P',p_ref,'T',T_ref,'Water') # Specific heat capacity of water at (p_ref, T_ref) [J/kg/K]
    dot_m_w = (1/(1+Xtot))*dot_m_v*Qc/(Cp_w*(T_cond_out-T_ref)) #[kg/s]: mass flow rate of condensed water

    print('Condensed water massflow: %.2f [kg/s]' %dot_m_w)

    h_cond_in = CP.PropsSI('H','P',p_ref,'T',T_ref,'Water')
    s_cond_in = CP.PropsSI('S','P',p_ref,'T',T_ref,'Water')
    e_cond_in = exergy(h_cond_in,s_cond_in)

    h_cond_out = CP.PropsSI('H','P',p_ref,'T',T_cond_out,'Water')
    s_cond_out = CP.PropsSI('S','P',p_ref,'T',T_cond_out,'Water')
    e_cond_out = exergy(h_cond_out,s_cond_out)

    eta_condex = dot_m_w*(e_cond_out-e_cond_in)/(((1+np.sum(X[:id_drum]))/(1+Xtot))*dot_m_v*eCond);

    eta_cyclen = eta_cyclen                                 #[-]: cycle energy efficiency     .489
    print('eta_cyclen: %.3f [-]' %eta_cyclen)
    eta_toten = P_e / (dot_m_f*LHV)                         #[-]: overall energy efficiency     .457
    print('eta_toten: %.3f [-]' %eta_toten)
    eta_cyclex = eta_cyclex                                 #[-]: cycle exergy efficiency     .849
    print('eta_cyclex: %.3f [-]' %eta_cyclex)
    eta_gen = eta_gen                                       #[-]: steam generator energy efficiency             .
    print('eta_gen: %.3f [-]' %eta_gen)
    eta_combex = dot_m_g*(e_f - e_r)/(dot_m_f*e_c)          #[-]: combustion exergy efficiency        .689
    print('eta_combex: %.3f [-]' %eta_combex)
    eta_chimnex = (e_f-e_exh) / (e_f-e_r)                   #[-]: chimney exergy efficiency     .991
    print('eta_chimnex: %.3f [-]' %eta_chimnex)
    eta_transex = dot_m_vG*( (1+Xtot)*(e_3-e_2) + (1+Xtot- X[-1])*(e_5-e_4) )/ (dot_m_g*(e_f-e_exh)) #[-]: bleedings heat exchangers overall exergy efficiency    .766
    print('eta_transex: %.3f [-]' %eta_transex)
    eta_gex = eta_transex*eta_chimnex*eta_combex            #[-]: steam generator exergy efficiency  .523
    print('eta_gex: %.3f [-]' %eta_gex)
    eta_totex = eta_gex*eta_cyclex*eta_mec                  #[-]: overall exergy efficiency  .440
    print('eta_totex: %.3f [-]' %eta_totex)
    eta_condex = eta_condex                                 #[-]: condenser exergy efficiency                           .
    print('eta_condex: %.3f [-]' %eta_condex)
    eta_rotex =  Wmtot/(emT-emP-emPe1-emPe2)    #[-]: pumps and turbines exergy efficiency  .918
    print('eta_rotex: %.3f [-]' %eta_rotex)
    print(75*'_')


    #
    #
    # ## Generate graph to export:
    # # ==========================
    #
    # # 1st figure : Energetic balance
    # fig_pie_en = plt.figure(1)
    # labels = 'Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Condensor loss \n'+'%.1f'%(loss_cond*1e-6)+' MW', 'Steam generator losses \n'+'%.1f'%(loss_gen*1e-6)+' MW'
    # sizes = [P_e*1e-6, loss_mec*1e-6, loss_cond*1e-6, loss_gen*1e-6]
    # plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=140)
    # plt.axis('equal')
    # plt.title("Primary energy flux " + "%.1f" %(LHV*dot_m_f*1e-6)+ " MW")
    #
    # # 2nd figure : Exergetic balance
    # fig_pie_ex = plt.figure(2)
    # labels = 'Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Condenser losses \n'+'%.1f'%(loss_condex*1e-6)+' MW', 'Turbine & pumps \n irreversibilities \n'+'%.1f'%(loss_rotex*1e-6)+' MW', 'Combustion \n irreversibilities \n'+'%.1f'%(loss_combex*1e-6)+' MW', 'Steam generator \n losses \n'+'%.1f'%((loss_gex-loss_combex-loss_chemex)*1e-6)+' MW', 'Chimney losses \n'+'%.1f'%(loss_chemex*1e-6)+' MW', 'Heat transfer irreversibilities \n in the feed-water heaters \n'+'%.1f'%(loss_transex*1e-6)+' MW'
    # sizes = [P_e*1e-6, loss_mec*1e-6, loss_condex*1e-6, loss_rotex*1e-6, loss_combex*1e-6, (loss_gex-loss_combex-loss_chemex)*1e-6, loss_chemex*1e-6, loss_transex*1e-6]
    # plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=30)
    # plt.axis('equal')
    # plt.title("Primary exergy flux " + "%.1f" %(e_c*dot_m_f*1e-6) + " MW")
    #
    # plt.show()

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
