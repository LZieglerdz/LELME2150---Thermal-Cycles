"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Final steam turbine function

@author: Laurent, Ziegler de Ziegleck aùf Rheingrüb, 03821500
         Dimitri, Boterberg                        , 10271700
@date: November 14, 2021
"""

import CoolProp
import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots.SimpleCycles import StateContainer
from CoolProp.Plots.Common import PropertyDict
import numpy as np
from thermochem import janaf
db = janaf.Janafdb();
import matplotlib.pyplot as plt

#
#===CHECK BEFORE SUBMISSION====================================================
#   - why is p_5 == 6200 [kPa] in the book (p91) and not 7000? -> 6200 may be hardcoded! -> don't forget to implement it correctly
#   - ditto for p_2 and p_3
#   - setting up the reference state; book seems to use the same one as the CoolProp library
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

def getCpBar(pi, pf, Ti, Tf, mix_comp, mix_conc):
    n=100
    cp = 0

    rangeT = np.linspace(Ti,Tf,n)
    rangeP = np.ones(n)*pi

    if Ti == Tf:
        return(cp)

    for i in range(n):
        cp += getCpMix(rangeT[i], rangeP[i], mix_comp, mix_conc) * (Tf-Ti) / (n*rangeT[i] * np.log(Tf/Ti))
    return(cp)


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
    # p_5 = 6200000
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
    # p_2 = 35000000
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

    dot_m_vc = P_e / ( eta_mec * Wmtot )
    dot_m_vG = dot_m_vc*(1+Xtot)

    dot_m_v = dot_m_vc                          #[kg/s]: mass flow rate of steam
    dot_m_f = dot_m_vG*Qh/(eta_gen*LHV)         #[kg/s]: mass flow rate of fuel
    dot_m_a = dot_m_f*lamb_ma1                  #[kg/s]: mass flow rate of air
    dot_m_g = dot_m_a+dot_m_f                   #[kg/s]: mass flow rate of flue gas

    # print('eps_exh: %.2f [-]' %eps_exh)
    # print('eps_p: %.2f [-]' %eps_p)
    # print('eta_gen: %.2f [-]' %eta_gen)
    # print('Qh: %.2f [J/kg]' %Qh)
    # print('LHV: %.2f [J/kg]' %LHV)
    # print('Steam massflow (cond.): %.2f [kg/s]' %dot_m_v)
    # print('Steam massflow (boil.): %.2f [kg/s]' %dot_m_vG)
    # print('Air massflow: %.2f [kg/s]' %dot_m_a)
    # print('Fuel massflow: %.2f [kg/s]' %dot_m_f)
    # print('Flue gas massflow: %.2f [kg/s]' %dot_m_g)
    # print(75*'_')

    #
    #===EFFICIENCIES===========================================================
    #
    e_r = 0
    c_pf = getMeanCp(p_ref, p_ref, T_ref, T_max_comb, R_f, comp, gas_prop)
    c_pf_int = getCpBar(p_ref, p_ref, T_ref, T_max_comb, comp, gas_prop)
    e_f = LHV/(lamb_ma1 + 1) - c_pf_int*T_ref*np.log(1 + LHV/((lamb_ma1 + 1)*c_pf*T_ref ) )
    e_exh = (getMeanCp(p_ref, p_ref, T_ref, T_max_comb, R_f, comp, gas_prop)*(T_exhaust-T_ref)-h_ref) - T_ref*(getMeanCp(p_ref, p_ref, T_ref, T_max_comb, R_f, comp, gas_prop)*np.log(T_exhaust/T_ref) - s_ref)


    # print(25*'_')
    # print('e_f (flue gasses after combustion): %.3f [kJ/kg]' %(e_f*1e-3))
    # print('e_exh (flue gasses at chimney): %.3f [kJ/kg]' %(e_exh*1e-3))
    # print(25*'_')


    Cp_w = CP.PropsSI('C','P',p_ref,'T',T_ref,'Water') # Specific heat capacity of water at (p_ref, T_ref) [J/kg/K]
    dot_m_w = (1/(1+Xtot))*dot_m_v*Qc/(Cp_w*(T_cond_out-T_ref)) #[kg/s]: mass flow rate of condensed water

    # print('Condensed water massflow: %.2f [kg/s]' %dot_m_w)

    h_cond_in = CP.PropsSI('H','P',p_ref,'T',T_ref,'Water')
    s_cond_in = CP.PropsSI('S','P',p_ref,'T',T_ref,'Water')
    e_cond_in = exergy(h_cond_in,s_cond_in)

    h_cond_out = CP.PropsSI('H','P',p_ref,'T',T_cond_out,'Water')
    s_cond_out = CP.PropsSI('S','P',p_ref,'T',T_cond_out,'Water')
    e_cond_out = exergy(h_cond_out,s_cond_out)

    eta_cyclen = eta_cyclen                                 #[-]: cycle energy efficiency     .489
    eta_toten = P_e / (dot_m_f*LHV)                         #[-]: overall energy efficiency     .457
    eta_gen = eta_gen                                       #[-]: steam generator energy efficiency             .

    eta_cyclex = eta_cyclex        #[-]: cycle exergy efficiency     .849
    eta_totex = P_e/ (dot_m_f*e_c)      #[-]: overall exergy efficiency  .440
    eta_transex = dot_m_vG*( (1+Xtot)*(e_3-e_2) + (1+Xtot- X[-1])*(e_5-e_4) )/ (dot_m_g*(e_f-e_exh)) #[-]: bleedings heat exchangers overall exergy efficiency    .766
    eta_combex = dot_m_g*(e_f - e_r)/(dot_m_f*e_c)      #[-]: combustion exergy efficiency        .689
    eta_gex = dot_m_vG*( (1+Xtot)*(e_3-e_2) + (1+Xtot- X[-1])*(e_5-e_4) )/(dot_m_f*e_c)     #[-]: steam generator exergy efficiency  .523
    eta_chimnex = eta_gex / (eta_transex * eta_combex)      #[-]: chimney exergy efficiency     .991
    eta_condex = dot_m_w*(e_cond_out-e_cond_in)/(((1+np.sum(X[:id_drum]))/(1+Xtot))*dot_m_v*eCond)  #[-]: condenser exergy efficiency                           .
    eta_rotex =  Wmtot/(emT-emP-emPe1-emPe2)                #[-]: pumps and turbines exergy efficiency  .918
    #
    print('EFFFICIENCIES')
    print('- energy')
    print('eta_cyclen: %.3f [-]  -- .489' %eta_cyclen)
    print('eta_toten: %.3f [-]   -- .457' %eta_toten)
    print('eta_gen: %.3f [-]     -- ' %eta_gen)
    print('- exergy')
    print('eta_cyclex: %.3f [-]  -- .849' %eta_cyclex)
    print('eta_combex: %.3f [-]  -- .689' %eta_combex)
    print('eta_chimnex: %.3f [-] -- .991' %eta_chimnex)
    print('eta_transex: %.3f [-] -- .766' %eta_transex)
    print('eta_gex: %.3f [-]     -- .523' %eta_gex)
    print('eta_totex: %.3f [-]   -- .440' %eta_totex)
    print('eta_condex: %.3f [-]  -- ' %eta_condex)
    print('eta_rotex: %.3f [-]   -- .918' %eta_rotex)
    print(75*'_')

    #
    #===LOSSES===========================================================
    #
    loss_gen = dot_m_f*LHV -  dot_m_vG*Qh
    loss_mec = dot_m_v*Wmtot - P_e
    loss_cond = (1 + np.sum(X[:id_drum])/(1 + Xtot))*dot_m_v*Qc

    loss_gex = dot_m_f*e_c - dot_m_vG*eSG
    loss_turbex = dot_m_v*(emT - WmT)
    loss_pumpex = dot_m_v*((WmP - emP) + (WmPe1 - emPe1) + (WmPe2 - emPe2))
    loss_rotex = loss_pumpex + loss_turbex
    loss_combex = dot_m_f*e_c - dot_m_g*e_f
    loss_condex = (1 + np.sum(X[:id_drum])/(1 + Xtot))*dot_m_v*eCond

    loss_transex = dot_m_g * (e_f-e_exh) - dot_m_vG* ( (1+Xtot)*(e_3-e_2) + (1+Xtot-X[-1])*(e_5-e_4) )
    loss_chemex = dot_m_g*(e_f-e_exh)

    loss_FWH = 0
    loss_FWH += (1+np.sum(X[:id_drum]))*((e7_b[1]-e_7)-(e9_b[0]-e_8));
    for i in range(id_drum-1):
        loss_FWH += X[i]*(e6_b[i+1]-e7_b[i+1])+np.sum(X[i:id_drum-1])*(e7_b[i+2]-e7_b[i+1])-(1+np.sum(X[:id_drum]))*(e9_b[i+1]-e9_b[i]);
    loss_FWH += X[id_drum-1]*(e6_b[id_drum]-e7_b[id_drum])-(1+np.sum(X[:id_drum]))*(e9_b[id_drum]-e9_b[id_drum-1]);
    for i in range(id_drum+1,nsout+reheat-1):
        loss_FWH += X[i]*(e6_b[i+1]-e7_b[i+1])+np.sum(X[i:nsout+reheat-1])*(e7_b[i+2]-e7_b[i+1])-(1+np.sum(X[id_drum+1:nsout+reheat]))*(e9_b[i+1]-e9_b[i]);
    loss_FWH += X[nsout+reheat-1]*(e6_b[nsout+reheat]-e7_b[nsout+reheat])-(1+np.sum(X[id_drum+1:nsout+reheat]))*(e9_b[nsout+reheat]-e9_b[nsout+reheat-1]);
    loss_FWH = loss_FWH*dot_m_v

    # "The bad way of computing things"
    loss_mec = P_e*(1-eta_mec)
    loss_chemex = (1-eta_chimnex) * dot_m_g * e_f
    loss_totex = (1-eta_totex) * dot_m_f * e_c
    # loss_transex = (1-eta_transex) * dot_m_g * (e_f-e_exh)
    # loss_FWH = loss_totex - (loss_gex+loss_condex+loss_rotex+loss_mec)


    print('LOSSES')
    print('loss_mec: %.3f [MW] -- Mechanical losses'      %(loss_mec*1e-6))
    print('loss_gen: %.3f [MW] -- Steam generator losses'      %(loss_gen*1e-6))
    print('loss_cond: %.3f [MW] -- Condensor loss'     %(loss_cond*1e-6))
    print(75*'_')
    print('loss_mec: %.3f [MW] -- 3MW -- Mechanical losses'      %(loss_mec*1e-6))
    print('loss_rotex: %.3f [MW] -- 26MW -- Turbine & pumps irreversibilities'    %(loss_rotex*1e-6))
    print('loss_gex: %.3f [MW] --   MW -- Steam generator losses -- combex+chemex+transex'      %(loss_gex*1e-6))
    print('loss_combex: %.3f [MW] -- 204MW -- Combustion irreversibilities'   %(loss_combex*1e-6))
    print('loss_chemex: %.3f [MW] -- 4MW -- Chimney losses'   %(loss_chemex*1e-6))
    print('loss_transex: %.3f [MW] -- 104MW -- Heat transfer irreversibilities in the steam generator'  %(loss_transex*1e-6))
    print('loss_FWH: %.3f [MW] -- 8MW -- Heat transfer irreversibilities in the feed-water heaters'    %(loss_FWH*1e-6))
    print('loss_condex: %.3f [MW] -- 18MW -- Condenser losses'   %(loss_condex*1e-6))
    print('loss_totex: %.3f [MW]'    %(loss_totex*1e-6))
    print(75*'_')



    ## Generate graph to export:
    # ==========================

    n = 50

    h_12 = np.linspace(h_1,h_2,n)
    s_12 = np.linspace(s_1,s_2,n)
    T_12 = np.zeros(n)
    for i in range(n):
        # s_12[i] = CP.PropsSI('S','P',p_2,'H',h_12[i],'Water')
        T_12[i] = CP.PropsSI('T','S',s_12[i],'H',h_12[i],'Water')

    T_23 = np.linspace(T_2,T_3,n)
    s_23 = np.zeros(n)
    h_23 = np.zeros(n)
    for i in range(n):
        s_23[i] = CP.PropsSI('S','P',p_3,'T',T_23[i],'Water')
        h_23[i] = CP.PropsSI('H','P',p_3,'T',T_23[i],'Water')

    p_34 = np.linspace(p_3,p_4,n)
    s_34 = np.zeros(n)
    h_34 = np.zeros(n)
    T_34 = np.zeros(n)
    for i in range(n):
        h_4is = CP.PropsSI('H','P',p_34[i],'S',s_3,'Water')
        h_34[i] = h_3 + (h_4is-h_3)*eta_is_HP
        s_34[i] = CP.PropsSI('S','P',p_34[i],'H',h_34[i],'Water')
        T_34[i] = CP.PropsSI('T','P',p_34[i],'H',h_34[i],'Water')

    T_45 = np.linspace(T_4,T_5,n)
    s_45 = np.zeros(n)
    h_45 = np.zeros(n)
    for i in range(n):
        s_45[i] = CP.PropsSI('S','P',p_4,'T',T_45[i],'Water')
        h_45[i] = CP.PropsSI('H','P',p_4,'T',T_45[i],'Water')

    p_56 = np.linspace(p_5,p_6,n)
    s_56 = np.zeros(n)
    h_56 = np.zeros(n)
    T_56 = np.zeros(n)
    for i in range(n):
        h_6is = CP.PropsSI('H','P',p_56[i],'S',s_5,'Water')
        h_56[i] = h_5 + (h_6is-h_5)*eta_is_LP
        s_56[i] = CP.PropsSI('S','P',p_56[i],'H',h_56[i],'Water')
        T_56[i] = CP.PropsSI('T','P',p_56[i],'H',h_56[i],'Water')

    T_67 = np.linspace(T_6,T_7,n)
    x_67 = np.linspace(x_6,x_7,n)
    s_67 = np.zeros(n)
    h_67 = np.zeros(n)
    for i in range(n):
        s_67[i] = CP.PropsSI('S','Q',x_67[i],'T',T_67[i],'Water')
        h_67[i] = CP.PropsSI('H','Q',x_67[i],'T',T_67[i],'Water')

    h_78 = np.linspace(h_7,h_8,n)
    s_78 = np.linspace(s_7,s_8,n)
    T_78 = np.zeros(n)
    for i in range(n):
        # s_12[i] = CP.PropsSI('S','P',p_2,'H',h_12[i],'Water')
        T_78[i] = CP.PropsSI('T','S',s_78[i],'H',h_78[i],'Water')

    T_89 = np.linspace(T_8,T_9,n)
    s_89 = np.zeros(n)
    h_89 = np.zeros(n)
    for i in range(n):
        s_89[i] = CP.PropsSI('S','P',p_9,'T',T_89[i],'Water')
        h_89[i] = CP.PropsSI('H','P',p_9,'T',T_89[i],'Water')

    h_97IV = np.linspace(h_9,h_7IV,n)
    s_97IV = np.linspace(s_9,s_7IV,n)
    T_97IV = np.zeros(n)
    for i in range(n):
        # s_12[i] = CP.PropsSI('S','P',p_2,'H',h_12[i],'Water')
        T_97IV[i] = CP.PropsSI('T','S',s_97IV[i],'H',h_97IV[i],'Water')

    h_7IV9IV = np.linspace(h_7IV,h_9IV,n)
    s_7IV9IV = np.linspace(s_7IV,s_9IV,n)
    T_7IV9IV = np.zeros(n)
    for i in range(n):
        # s_12[i] = CP.PropsSI('S','P',p_2,'H',h_12[i],'Water')
        T_7IV9IV[i] = CP.PropsSI('T','S',s_7IV9IV[i],'H',h_7IV9IV[i],'Water')

    h_9IV9VIII = np.linspace(h_9IV,h_9VIII,n)
    s_9IV9VIII = np.linspace(s_9IV,s_9VIII,n)
    T_9IV9VIII = np.zeros(n)
    for i in range(n):
        # s_12[i] = CP.PropsSI('S','P',p_2,'H',h_12[i],'Water')
        T_9IV9VIII[i] = CP.PropsSI('T','S',s_9IV9VIII[i],'H',h_9IV9VIII[i],'Water')

    h_9VIII1 = np.linspace(h_9VIII,h_1,n)
    s_9VIII1 = np.linspace(s_9VIII,s_1,n)
    T_9VIII1 = np.zeros(n)
    for i in range(n):
        # s_12[i] = CP.PropsSI('S','P',p_2,'H',h_12[i],'Water')
        T_9VIII1[i] = CP.PropsSI('T','S',s_9VIII1[i],'H',h_9VIII1[i],'Water')


    # BLEEDINGS
    h_6I7I = np.linspace(h_6I,h_7I,n, endpoint=False)
    s_6I7I = np.zeros(n)
    T_6I7I = np.zeros(n)
    for i in range(n):
        s_6I7I[i] = CP.PropsSI('S','P',p_6I,'H',h_6I7I[i],'Water')
        T_6I7I[i] = CP.PropsSI('T','P',p_6I,'H',h_6I7I[i],'Water')

    h_6II7II = np.linspace(h_6II,h_7II,n, endpoint=False)
    s_6II7II = np.zeros(n)
    T_6II7II = np.zeros(n)
    for i in range(n):
        s_6II7II[i] = CP.PropsSI('S','P',p_6II,'H',h_6II7II[i],'Water')
        T_6II7II[i] = CP.PropsSI('T','P',p_6II,'H',h_6II7II[i],'Water')

    h_6III7III = np.linspace(h_6III,h_7III,n, endpoint=False)
    s_6III7III = np.zeros(n)
    T_6III7III = np.zeros(n)
    for i in range(n):
        s_6III7III[i] = CP.PropsSI('S','P',p_6III,'H',h_6III7III[i],'Water')
        T_6III7III[i] = CP.PropsSI('T','P',p_6III,'H',h_6III7III[i],'Water')

    h_6IV7IV = np.linspace(h_6IV,h_7IV,n, endpoint=False)
    s_6IV7IV = np.zeros(n)
    T_6IV7IV = np.zeros(n)
    for i in range(n):
        s_6IV7IV[i] = CP.PropsSI('S','P',p_6IV,'H',h_6IV7IV[i],'Water')
        T_6IV7IV[i] = CP.PropsSI('T','P',p_6IV,'H',h_6IV7IV[i],'Water')

    h_6V7V = np.linspace(h_6V,h_7V,n, endpoint=False)
    s_6V7V = np.zeros(n)
    T_6V7V = np.zeros(n)
    for i in range(n):
        s_6V7V[i] = CP.PropsSI('S','P',p_6V,'H',h_6V7V[i],'Water')
        T_6V7V[i] = CP.PropsSI('T','P',p_6V,'H',h_6V7V[i],'Water')

    h_6VI7VI = np.linspace(h_6VI,h_7VI,n, endpoint=False)
    s_6VI7VI = np.zeros(n)
    T_6VI7VI = np.zeros(n)
    for i in range(n):
        s_6VI7VI[i] = CP.PropsSI('S','P',p_6VI,'H',h_6VI7VI[i],'Water')
        T_6VI7VI[i] = CP.PropsSI('T','P',p_6VI,'H',h_6VI7VI[i],'Water')

    h_6VII7VII = np.linspace(h_6VII,h_7VII,n, endpoint=False)
    s_6VII7VII = np.zeros(n)
    T_6VII7VII = np.zeros(n)
    for i in range(n):
        s_6VII7VII[i] = CP.PropsSI('S','P',p_6VII,'H',h_6VII7VII[i],'Water')
        T_6VII7VII[i] = CP.PropsSI('T','P',p_6VII,'H',h_6VII7VII[i],'Water')

    h_6VIII7VIII = np.linspace(h_6VIII,h_7VIII,n, endpoint=False)
    s_6VIII7VIII = np.zeros(n)
    T_6VIII7VIII = np.zeros(n)
    for i in range(n):
        s_6VIII7VIII[i] = CP.PropsSI('S','P',p_6VIII,'H',h_6VIII7VIII[i],'Water')
        T_6VIII7VIII[i] = CP.PropsSI('T','P',p_6VIII,'H',h_6VIII7VIII[i],'Water')

    s_7I7 = np.ones(n)*s_7I
    p_7I7 = np.linspace(p_7I,p_7,n)
    h_7I7 = np.zeros(n)
    T_7I7 = np.zeros(n)
    for i in range(n):
        h_7I7[i] = CP.PropsSI('H','S',s_7I7[i],'P',p_7I7[i],'Water')
        T_7I7[i] = CP.PropsSI('T','S',s_7I7[i],'P',p_7I7[i],'Water')

    s_7II7I = np.ones(n)*s_7II
    p_7II7I = np.linspace(p_7II,p_7I,n)
    h_7II7I = np.zeros(n)
    T_7II7I = np.zeros(n)
    for i in range(n):
        h_7II7I[i] = CP.PropsSI('H','S',s_7II7I[i],'P',p_7II7I[i],'Water')
        T_7II7I[i] = CP.PropsSI('T','S',s_7II7I[i],'P',p_7II7I[i],'Water')

    s_7III7II = np.ones(n)*s_7III
    p_7III7II = np.linspace(p_7III,p_7II,n)
    h_7III7II = np.zeros(n)
    T_7III7II = np.zeros(n)
    for i in range(n):
        h_7III7II[i] = CP.PropsSI('H','S',s_7III7II[i],'P',p_7III7II[i],'Water')
        T_7III7II[i] = CP.PropsSI('T','S',s_7III7II[i],'P',p_7III7II[i],'Water')

    s_7IV7III = np.ones(n)*s_7IV
    p_7IV7III = np.linspace(p_7IV,p_7III,n)
    h_7IV7III = np.zeros(n)
    T_7IV7III = np.zeros(n)
    for i in range(n):
        h_7IV7III[i] = CP.PropsSI('H','S',s_7IV7III[i],'P',p_7IV7III[i],'Water')
        T_7IV7III[i] = CP.PropsSI('T','S',s_7IV7III[i],'P',p_7IV7III[i],'Water')

    s_7V7IV = np.ones(n)*s_7V
    p_7V7IV = np.linspace(p_7V,p_7IV,n)
    h_7V7IV = np.zeros(n)
    T_7V7IV = np.zeros(n)
    for i in range(n):
        h_7V7IV[i] = CP.PropsSI('H','S',s_7V7IV[i],'P',p_7V7IV[i],'Water')
        T_7V7IV[i] = CP.PropsSI('T','S',s_7V7IV[i],'P',p_7V7IV[i],'Water')

    s_7VI7V = np.ones(n)*s_7VI
    p_7VI7V = np.linspace(p_7VI,p_7V,n)
    h_7VI7V = np.zeros(n)
    T_7VI7V = np.zeros(n)
    for i in range(n):
        h_7VI7V[i] = CP.PropsSI('H','S',s_7VI7V[i],'P',p_7VI7V[i],'Water')
        T_7VI7V[i] = CP.PropsSI('T','S',s_7VI7V[i],'P',p_7VI7V[i],'Water')

    s_7VII7VI = np.ones(n)*s_7VII
    p_7VII7VI = np.linspace(p_7VII,p_7VI,n)
    h_7VII7VI = np.zeros(n)
    T_7VII7VI = np.zeros(n)
    for i in range(n):
        h_7VII7VI[i] = CP.PropsSI('H','S',s_7VII7VI[i],'P',p_7VII7VI[i],'Water')
        T_7VII7VI[i] = CP.PropsSI('T','S',s_7VII7VI[i],'P',p_7VII7VI[i],'Water')

    s_7VIII7VII = np.ones(n)*s_7VIII
    p_7VIII7VII = np.linspace(p_7VIII,p_7VII,n)
    h_7VIII7VII = np.zeros(n)
    T_7VIII7VII = np.zeros(n)
    for i in range(n):
        h_7VIII7VII[i] = CP.PropsSI('H','S',s_7VIII7VII[i],'P',p_7VIII7VII[i],'Water')
        T_7VIII7VII[i] = CP.PropsSI('T','S',s_7VIII7VII[i],'P',p_7VIII7VII[i],'Water')

    s_points = np.append(s_12, [s_23,s_34,s_45,s_56,s_67,s_78,s_89,s_97IV,s_7IV9IV,s_9IV9VIII] )#*1e-3
    h_points = np.append(h_12, [h_23,h_34,h_45,h_56,h_67,h_78,h_89,h_97IV,h_7IV9IV,h_9IV9VIII] )#*1e-3
    T_points = np.append(T_12, [T_23,T_34,T_45,T_56,T_67,T_78,T_89,T_97IV,T_7IV9IV,T_9IV9VIII] )
    T_points -= 273.15

    s_bI = np.append(s_6I7I ,s_7I7 )#*1e-
    s_bII = np.append(s_6II7II,s_7II7I )
    s_bIII = np.append(s_6III7III,s_7III7II )
    s_bIV = np.append(s_6IV7IV,s_7IV7III )
    s_bV = np.append(s_6V7V,s_7V7IV )
    s_bVI = np.append(s_6VI7VI,s_7VI7V )
    s_bVII = np.append(s_6VII7VII,s_7VII7VI )
    s_bVIII = np.append(s_6VIII7VIII,s_7VIII7VII )

    h_bI = np.append(h_6I7I, h_7I7 )#*1e-
    h_bII = np.append(h_6II7II,h_7II7I )
    h_bIII = np.append(h_6III7III,h_7III7II )
    h_bIV = np.append(h_6IV7IV,h_7IV7III )
    h_bV = np.append(h_6V7V,h_7V7IV )
    h_bVI = np.append(h_6VI7VI,h_7VI7V )
    h_bVII = np.append(h_6VII7VII,h_7VII7VI )
    h_bVIII = np.append(h_6VIII7VIII,h_7VIII7VII )

    T_bI = np.append(T_6I7I, T_7I7 )#*1e-
    T_bII = np.append(T_6II7II,T_7II7I )
    T_bIII = np.append(T_6III7III,T_7III7II )
    T_bIV = np.append(T_6IV7IV,T_7IV7III )
    T_bV = np.append(T_6V7V,T_7V7IV )
    T_bVI = np.append(T_6VI7VI,T_7VI7V )
    T_bVII = np.append(T_6VII7VII,T_7VII7VI )
    T_bVIII = np.append(T_6VIII7VIII,T_7VIII7VII )

    T_bI -= 273.15
    T_bII -= 273.15
    T_bIII -= 273.15
    T_bIV -= 273.15
    T_bV -= 273.15
    T_bVI -= 273.15
    T_bVII -= 273.15
    T_bVIII -= 273.15

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

    # 1st figure : Energetic balance
    fig_pie_en = plt.figure(1)
    labels = 'Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Condensor loss \n'+'%.1f'%(loss_cond*1e-6)+' MW', 'Steam generator losses \n'+'%.1f'%(loss_gen*1e-6)+' MW'
    sizes = [P_e*1e-6, loss_mec*1e-6, loss_cond*1e-6, loss_gen*1e-6]
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')
    plt.title("Primary energy flux " + "%.1f" %(LHV*dot_m_f*1e-6)+ " MW")


    # 2nd figure : Exergetic balance
    fig_pie_ex = plt.figure(2)
    labels = ['Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW',    'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW',    'Condenser losses\n'+'%.1f'%(loss_condex*1e-6)+' MW',    'Turbine & pumps \n irreversibilities \n'+'%.1f'%(loss_rotex*1e-6)+' MW', 'Heat transfer irreversibilities \n in the steam generator \n'+'%.1f'%(loss_transex*1e-6)+' MW','Heat transfer irreversibilities \n in the feed-water heaters \n'+'%.1f'%(loss_FWH*1e-6)+' MW',    'Combustion \n irreversibilities \n'+'%.1f'%(loss_combex*1e-6)+' MW','Chimney losses \n'+'%.1f'%(loss_chemex*1e-6)+' MW']
    sizes = [P_e*1e-6, loss_mec*1e-6, loss_condex*1e-6, loss_rotex*1e-6, loss_transex*1e-6, loss_FWH*1e-6,loss_combex*1e-6,loss_chemex*1e-6]
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=30)
    plt.axis('equal')
    plt.title("Primary exergy flux " + "%.1f" %(e_c*dot_m_f*1e-6) + " MW")


    # 3rd figure: T-s
    fig_Ts_diagram = plt.figure(3)
    fig_Ts_diagram = PropertyPlot('water', 'Ts', unit_system='EUR')
    fig_Ts_diagram.calc_isolines(CoolProp.iQ, num=11)
    fig_Ts_diagram.set_axis_limits([0., 9, 0, T_max+100-273.15])
    plt.grid(True)
    plt.plot(s_points*1e-3, T_points, c="red")
    plt.plot(s_bI*1e-3, T_bI, c="red")
    plt.plot(s_bII*1e-3, T_bII, c="red")
    plt.plot(s_bIII*1e-3, T_bIII, c="red")
    plt.plot(s_bIV*1e-3, T_bIV, c="red")
    plt.plot(s_bV*1e-3, T_bV, c="red")
    plt.plot(s_bVI*1e-3, T_bVI, c="red")
    plt.plot(s_bVII*1e-3, T_bVII, c="red")
    plt.plot(s_bVIII*1e-3, T_bVIII, c="red")
    plt.plot(np.array(s)*1e-3, np.array(T)-273.15, 'ko')
    plt.title('T-s diagram of the cycle')
    plt.xlabel("s $[kJ/kg/K]$")
    plt.ylabel("T $[°C]$")


    # 4th figure: h-s
    fig_hs_diagram = plt.figure(4)
    fig_hs_diagram = PropertyPlot('water', 'HS', unit_system='EUR')
    fig_hs_diagram.calc_isolines(CoolProp.iQ, num=11)
    fig_hs_diagram.set_axis_limits([0., 9, 0, 4000])
    plt.grid(True)
    plt.plot(s_points*1e-3, h_points*1e-3, c="red")
    plt.plot(s_bI*1e-3, h_bI*1e-3, c="red")
    plt.plot(s_bII*1e-3, h_bII*1e-3, c="red")
    plt.plot(s_bIII*1e-3, h_bIII*1e-3, c="red")
    plt.plot(s_bIV*1e-3, h_bIV*1e-3, c="red")
    plt.plot(s_bV*1e-3, h_bV*1e-3, c="red")
    plt.plot(s_bVI*1e-3, h_bVI*1e-3, c="red")
    plt.plot(s_bVII*1e-3, h_bVII*1e-3, c="red")
    plt.plot(s_bVIII*1e-3, h_bVIII*1e-3, c="red")
    plt.plot(np.array(s)*1e-3, np.array(h)*1e-3, 'ko')
    plt.title('h-s diagram of the cycle')
    plt.xlabel("s $[kJ/kg/K]$")
    plt.ylabel("h $[kJ/kg]$")

    fig_pie_en.savefig('pie_en_diag.png', dpi=200)
    fig_pie_ex.savefig('pie_ex_diag.png', dpi=200)
    fig_Ts_diagram.savefig('Ts_diag.png', dpi=200)
    fig_hs_diagram.savefig('hs_diag.png', dpi=200)

    if display:
        plt.show()

    DAT = (p,T,h,s,e,x)
    COMBUSTION = (LHV,e_c,excess_air,cp_gas,gas_prop)
    MASSFLOW = (dot_m_a,dot_m_f,dot_m_g,dot_m_v)
    ETA = (eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex)
    DATEN = (loss_mec,loss_gen,loss_cond)
    DATEX = (loss_mec,loss_rotex,loss_combex,loss_chemex,loss_transex,loss_totex,loss_condex)
    XMASSFLOW = (x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII)
    FIG = (fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram)
    out = (ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG)
    return out
