"""
LELME2150 - Thermal cycles
Homework 4

GST function

@author: Laurent, Ziegler de Ziegleck aùf Rheingrüb, 03821500
         Dimitri, Boterberg                        , 10271700
@date: December 15, 2021
"""

#
#===IMPORT PACKAGES============================================================
#

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
#===GLOBAL VARIABLES===========================================================
#
T_S = 273.153         # [K]
p_S = 1e5             # [Pa]

# air composition
N2_conc = .7808         # %mol
O2_conc = .2095         # %mol
H2O_conc = 0          # %mol
CO2_conc = .0004          # %mol

Mm_N2 = CP.PropsSI("MOLARMASS", "N2")        # kg/mol
Mm_O2 = CP.PropsSI("MOLARMASS", "O2")        # kg/mol
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

def get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3,iter, lam_est):
    #  Combustion supposée complète
    w = lam_est*(z + (y-2*x)/4)
    a_0 = (x + 2*w - (y/2) - 2*z)/2
    Mm_f, flue_conc_mol, R_f = FlueGasFrac(x, y, lam_est)
    flue_conc_mass = np.multiply(flue_conc_mol, [Mm_N2, Mm_O2, Mm_CO2, Mm_H2O])/Mm_f

    cp_3S = getMeanCp(p_S, p_3, T_S, T_3, R_f, comp, flue_conc_mass)
    cp_32 = getMeanCp(p_2, p_3, T_2, T_3, R_f, comp, flue_conc_mass)

    lam = (get_LHV(x,y)-cp_3S*(T_3-T_S))/ (get_ma1(x,y)*cp_32*(T_3-T_2))

    if iter == 0:
        print('The function get_lambda did not converge.')
        return(-1)
    if np.abs(lam-lam_est) <= 1e-6 :
        return(lam)
    else:
        return( get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3, iter-1, lam) )

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
    x_f = [x_N2f,x_O2f,x_CO2f,x_H2Of];
    x_O2f = x_O2f/sum(x_f);
    x_N2f = x_N2f/sum(x_f);
    x_CO2f = x_CO2f/sum(x_f);
    x_H2Of = x_H2Of/sum(x_f);
    Mm_f = Mm_N2*x_N2f + Mm_O2*x_O2f + Mm_CO2*x_CO2f + Mm_H2O*x_H2Of; # Molar mass of fluegas [g/mol]
    R_f = janaf.R/Mm_f; # Ideal gas constant (R*) for exhaust gas [J/kg/K]

    return [Mm_f, [x_N2f, x_O2f, x_CO2f, x_H2Of], R_f]

def getCpMix(T, p, mix_comp, mix_conc):
    sum = 0
    for i in range(len(mix_comp)):
         sum += (mix_conc[i]*CP.PropsSI('CPMASS', 'T', T, 'P', p, mix_comp[i]) )
    return(sum)

def getMeanCp(p_in, p_out, T_in, T_out, R, mix_comp, mix_conc):               # Gives the mean specific heat capacity
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
#===POLYTROPIC TEMPERATURE=====================================================
#
def getPolytropicTemp(p_in, p_out, T_in, T_out, R, eta_pi, iter, mix_comp, mix_conc):
    if iter < 0:
        print("Function does not converge.")
        return(-1)
    cp = getMeanCp(p_in,p_out,T_in,T_out, R, mix_comp, mix_conc )
    T = T_in*(p_out/p_in)**(R/(eta_pi*cp))
    # print(iter, T-273.15, cp)
    if np.abs(T_out-T) < 1e-12:
        return(T)
    else:
        return(getPolytropicTemp(p_in, p_out, T_in, T, R, eta_pi, iter-1, mix_comp, mix_conc))

#
#===GST CYCLE================================================================
#

# SPECIFICATIONS:
#   GST(P_e,options,display) computes the thermodynamics states for
#   the combined power plant of pp. 166, fig 4.19
#   based on several inputs (given in options) and based on a given electricity
#   production P_e.  It returns the main results. It can also plot the T-s and
#   h-s as well as the energy and exergy pies if display == True
# INPUTS (/!\ some inputs can be dependent on others /!\):
#   - P_eg              [W]: Net power output of gas turbine
#   - P_es              [W]: Net power output of steam turbine
#   - options: tuple containing the parametric values of the GT & ST
#       o p_1g          [Pa]: inlet pressure (ambient)
#       o T_1g          [K]: inlet temperature (ambient)
#       o r             [-]: compression ratio
#       o T_3g          [K]: combustor outlet temperature
#       o k_cc          [-]: pressure losses coeffcient
#       o eta_pi_c      [-]: compressor polytropic efficiency
#       o eta_pi_t      [-]: turbine polytropic efficiency
#       o k_mec         [-]: shaft losses
#       o comb          [na]: tuple containing combustion data
#           * x         [-]: O_x/C ratio of fuel (e.g. 0.0 in CH_4)
#           * y         [-]: H_y/C ratio of fuel (e.g. 4.0 in CH_4)
#
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
#       o eta_mec       [-]: shafts bearings mechanical efficiency
#       o eta_pump      [-]: internal efficiency of the pump (see text book pp. 53)
#       o eta_is_turb   [na]: tuple containing turbines isentropic efficiencies (see text book pp. 55)
#           * eta_is_HP [-]: high pressure turbine isentropic efficiency
#           * eta_is_LP [-]: low pressure turbine isentropic efficiency
#   - display: bool to choose to plot the T-s & h-s diagrams and the energy and exergy pie charts (True or False)
#
# OUTPUTS: tuple containing...

#   - DAT_g: tuple containing the GT cycle state data
#      o p_g [Pa]: tuple containing the pressure at each state
#      o T_g [T]: tuple containing the temperature at each state
#      o h_g [J/kg]: tuple containing the enthalpy at each state
#      o s_g [J/kg/K]: tuple containing the entropy at each state
#      o e_g [J/kg]: tuple containing the exergy at each state
#   - COMBUSTION: tuple containing the combustion parameters
#      o LHV [J/kg]: the fuel Lower Heating Value
#      o e_c [J/kg]: the fuel exergy
#      o excess_air [-]: the excess air of the combustion
#      o cp_gas [J/kg/K]: the flue gas specific heat capacity at the combustor outlet
#      o gas_prop: list containing the proportion of ['CO2','H2O','N2','O2'] in the flue gas respectively
#   - MASSFLOW_g: tuple containing the massflow rates
#      o dot_m_a [kg/s]: mass flow rate of air
#      o dot_m_f [kg/s]: mass flow rate of fuel
#      o dot_m_g [kg/s]: mass flow rate of flue gas
#   - ETA_g: tuple containing the efficiencies (see text book pp. 113-125)
#       o eta_cyclen [-]: cycle energy efficiency
#       o eta_toten [-]: overall energy efficiency
#       o eta_cyclex [-]: cycle exergy efficiency
#       o eta_totex [-]: overall exergy efficiency
#       o eta_rotex [-]: compressor-turbine exergy efficiency
#       o eta_combex [-]: combustion exergy efficiency

#   - DAT: tuple containing the ST cycle state data (1,2,3,4,5,6,7,8,8',8",9,9',9",10,10',10")
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
#   - FIG: tuple containing the figures to be diplayed
#       o fig_pie_en: pie chart of energy losses
#       o fig_pie_ex: pie chart of exergy losses
#       o fig_Ts_diagram: T-s diagram of the GT cycle
#       o fig_hs_diagram: h-s diagram of the GT cycle
def GST(P_eg, P_es, options, display):

    (p_1g,p_2g,p_3g,p_4g,p_5g,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_8p,p_8pp,p_9,p_9p,p_9pp,p_10p,p_10pp) = np.zeros(20)
    (T_1g,T_2g,T_3g,T_4g,T_5g,T_1,T_2,T_3,T_4,T_5,T_6,T_7,T_8,T_8p,T_8pp,T_9,T_9p,T_9pp,T_10p,T_10pp) = np.zeros(20)
    (s_1g,s_2g,s_3g,s_4g,s_5g,s_1,s_2,s_3,s_4,s_5,s_6,s_7,s_8,s_8p,s_8pp,s_9,s_9p,s_9pp,s_10p,s_10pp) = np.zeros(20)
    (h_1g,h_2g,h_3g,h_4g,h_5g,h_1,h_2,h_3,h_4,h_5,h_6,h_7,h_8,h_8p,h_8pp,h_9,h_9p,h_9pp,h_10p,h_10pp) = np.zeros(20)
    (e_1g,e_2g,e_3g,e_4g,e_5g,e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_8p,e_8pp,e_9,e_9p,e_9pp,e_10p,e_10pp) = np.zeros(20)
    (x_1g,x_2g,x_3g,x_4g,x_5g,x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_8p,x_8pp,x_9,x_9p,x_9pp,x_10p,x_10pp) = np.zeros(20)

    # Process input variables (p.167 book)--------------------------------------------------
    T_ref = 288.15;
    p_ref = 1e+5;
    T_1g = 15+273.15;   # Gas turbine air intake temperature [K]
    p_1g = 100e3;       # [Pa]
    T_3g = 1250+273.15; # Gas turbine flue gas temperature [K]
    r = 15;             # Compression ratio [-]
    eta_pi_c = 0.9;
    eta_pi_t = 0.9;
    k_cc = 1-0.05;
    k_mec = 0.015;

    x = 0; y = 4; # CH4

    T_pinch_approach = 50;        # Approach pinch temperature [K]
    T_pinch = 10;                 # Pinch temperature [K]
    eta_pump = 0.85;              # Extraction pump (1->2) efficiency [-]
    eta_is_HP = 0.92;             # HP turbine isentropic efficiency [-]
    eta_is_LP = 0.90;             # LP turbine isentropic efficiency [-]
    p_3 = 11e6;                   # HP inlet pressure [Pa] (Maximum steam pressure)
    p_5 = 2.8e6;                  # IP inlet pressure [Pa] (Reaheating pressure)
    T_5 = 565+273.15;             # IP inlet temperature [K] (Reheating temperature)
    p_6 = 0.4e6;                  # LP inlet pressure [Pa]
    T_6 = 318+273.15;             # LP inlet temperature [K]
    p_7 = 6e3;                    # LP outlet pressure [Pa]
    x_7 = 0.95;                   # LP outlet steam quality (condenser inlet) = minimum possible vapor quality after final expansion [-]
    # T_1 = T_cond_out ?
    # T_pinch_cond ?

    #===Ref. State for steam turbine============================================

    T_ref = T_ref   #[K]
    p_ref = p_ref   #[Pa]

    # Dmolar = CP.PropsSI("Dmolar", "T", 273.15+36.16, "P", 6000, 'Water')
    # CP.set_reference_state('Water', T_ref, Dmolar, 0, 0)

    # comp = ["N2","O2","CO2","H2O"]
    # air_conc = [N2_conc,O2_conc,CO2_conc,H2O_conc]


    h_ref = CP.PropsSI('H','P',p_ref,'T', T_ref,'Water') #[J/kg]
    s_ref = CP.PropsSI('S','P',p_ref,'T',T_ref,'Water') #[J/kgK]
    h_setref = 15.2e3
    s_setref = .054e3

    for i in np.arange(len(comp[:-1])):
        elem = comp[i]
        Dmolar = CP.PropsSI("Dmolar", "T", 273.15, "P", p_1g, elem)
        CP.set_reference_state(elem, T_ref-15, Dmolar, 0, 0)


    #===Gas Cycle=============================================================
    ma1 = get_ma1(x,y)
    LHV = get_LHV(x, y)                         #[J/kg]: the fuel Lower Heating Value                                  
    cp_gas, e_c  = CombEx(x, y)                 #[J/kg/K]: the flue gas specific heat capacity at the combustor outlet, [J/kg]: the fuel exergy
    f = e_c/LHV


    R_Star = 0
    for i in range(len(comp)):
        R_Star += air_conc[i]*CP.PropsSI('GAS_CONSTANT','T', T_1g,'P', p_1g, comp[i])/Mm[i] # [J/kgK]

    # State 1g
    s_1g = .79*CP.PropsSI('S','T', T_1g,'P', p_1g, "N2") + .21*CP.PropsSI('S','T', T_1g,'P', p_1g, "O2")
    h_1g = .79*CP.PropsSI('H','T', T_1g,'P', p_1g, "N2") + .21*CP.PropsSI('H','T', T_1g,'P', p_1g, "O2")
    e_1g = (h_1g - h_1g) - T_1g*(s_1g - s_1g)

    # State 2g -- 1g->2g: Polytropic Compression
    p_2g = r*p_1g
    T_2g = getPolytropicTemp(p_1g, p_2g, T_1g, T_1g, R_Star, eta_pi_c, 100, comp, air_conc)
    cp_2g = getMeanCp(p_1g,p_2g,T_1g,T_2g, R_Star,comp, air_conc)
    h_2g = h_1g + cp_2g*(T_2g-T_1g)
    s_2g = s_1g + cp_2g*np.log(T_2g/T_1g) - R_Star*np.log(p_2g/p_1g)
    e_2g = (h_2g - h_1g) - T_1g*(s_2g - s_1g)

    # State 3g -- 2g->3g: Isobar Combustion
    p_3g = k_cc*p_2g
    lamb = get_lambda('CH4', 1, y, x, T_2g, T_3g, p_2g, p_3g, 200, 1)
    excess_air = lamb
    print("lambda: %.2f" %lamb)

    Mm_f, flue_conc_mol, R_f = FlueGasFrac(x, y, lamb)
    flue_conc_mass = np.multiply(flue_conc_mol, [Mm_O2, Mm_N2, Mm_CO2, Mm_H2O])/Mm_f
    cp_3g = getMeanCp(p_2g,p_3g,T_2g,T_3g, R_f, comp, flue_conc_mass)
    h_3g = h_2g + cp_3g*(T_3g-T_2g)
    s_3g = s_2g + cp_3g*np.log(T_3g/T_2g) - R_f*np.log(p_3g/p_2g)
    e_3g = (h_3g - h_1g) - T_1g*(s_3g - s_1g)

    lamb_ma1 = lamb*ma1;

    # State 4g -- 3g->4g: Polytropic Expansion
    p_4g = p_1g
    T_4g = getPolytropicTemp(p_3g, p_4g, T_3g, T_3g, R_f , 1/eta_pi_t, 100, ['N2','O2','CO2','H2O'], flue_conc_mass)
    cp_4g = getMeanCp(p_3g,p_4g,T_4g,T_3g, R_f, comp, flue_conc_mass)
    h_4g = h_3g + cp_4g*(T_4g-T_3g)
    s_4g = s_3g + cp_4g*np.log(T_4g/T_3g) - R_f*np.log(p_4g/p_3g)
    e_4g = (h_4g - h_1g) - T_1g*(s_4g - s_1g)


    eta_mecg = 1 - k_mec*((1+(1/lamb_ma1))*(h_3g - h_4g) + (h_2g - h_1g))/((1+(1/lamb_ma1))*(h_3g - h_4g) - (h_2g - h_1g)); #Mechanical efficiency


    ## Gas cycle Mass flows
    # =====================

    dotm_a = P_eg/(eta_mecg*((1+(1/lamb_ma1))*(h_3g - h_4g) - (h_2g - h_1g))); # Air mass flow rate [kg_air/s]
    dotm_f = dotm_a/lamb_ma1;                                                  # Combustible mass flow rate [kg_comb/s]
    dotm_g = (1+(1/lamb_ma1))*dotm_a;                                          # Fluegas mass flow rate [kg_fluegas/s]

    gas_prop = [flue_conc_mass[0]*(1+(1/lamb_ma1))*dotm_a, # N2
                flue_conc_mass[1]*(1+(1/lamb_ma1))*dotm_a, # O2
                flue_conc_mass[2]*(1+(1/lamb_ma1))*dotm_a, # CO2
                flue_conc_mass[3]*(1+(1/lamb_ma1))*dotm_a] # H2O

    #===Steam Cycle=============================================================


    def exergy(h,s):
        return( (h-h_ref) - T_ref*(s-s_ref) )


    # State 3 -> HP inlet:
    p_3 = p_3
    T_3 = T_4g - T_pinch_approach
    h_3 = CP.PropsSI('H','P',p_3,'T',T_3,'Water')
    s_3 = CP.PropsSI('S','P',p_3,'T',T_3,'Water')
    x_3 = CP.PropsSI('Q','P',p_3,'T',T_3,'Water')
    e_3 = exergy(h_3,s_3)

    # State 10' -> Inlet Evap HP (10'->10"->3 isobar):
    p_10p = p_3
    x_10p = 0 # Saturated liquid
    T_10p = CP.PropsSI('T','P',p_10p,'Q',x_10p,'Water') # = T_6 = T_9
    h_10p = CP.PropsSI('H','P',p_10p,'Q',x_10p,'Water')
    s_10p = CP.PropsSI('S','P',p_10p,'Q',x_10p,'Water')
    e_10p = exergy(h_10p,s_10p)

    # State 10" -> Outlet Evap HP (10'->10"->3 isobar):
    p_10pp = p_10p
    x_10pp = 1 # Saturated vapor
    T_10pp = CP.PropsSI('T','P',p_10pp,'Q',x_10pp,'Water') # = T_6 = T_9
    h_10pp = CP.PropsSI('H','P',p_10pp,'Q',x_10pp,'Water')
    s_10pp = CP.PropsSI('S','P',p_10pp,'Q',x_10pp,'Water')
    e_10pp = exergy(h_10pp,s_10pp)

    # State 5 -> IP inlet:
    p_5 = p_5
    T_5 = T_5
    h_5 = CP.PropsSI('H','P',p_5,'T',T_5,'Water')
    s_5 = CP.PropsSI('S','P',p_5,'T',T_5,'Water')
    x_5 = CP.PropsSI('Q','P',p_5,'T',T_5,'Water')
    e_5 = exergy(h_5,s_5)

    # State 9' -> Inlet Evap IP (4->5->9->9'->9" isobar):
    p_9p = p_5
    x_9p = 0 # Saturated liquid
    T_9p = CP.PropsSI('T','P',p_9p,'Q',x_9p,'Water')
    h_9p = CP.PropsSI('H','P',p_9p,'Q',x_9p,'Water')
    s_9p = CP.PropsSI('S','P',p_9p,'Q',x_9p,'Water')
    e_9p = exergy(h_9p,s_9p)

    # State 9" -> Outlet Evap IP (4->5->9->9'->9" isobar):
    p_9pp = p_9p
    x_9pp = 1 # Saturated vapor
    T_9pp = CP.PropsSI('T','P',p_9pp,'Q',x_9pp,'Water')
    h_9pp = CP.PropsSI('H','P',p_9pp,'Q',x_9pp,'Water')
    s_9pp = CP.PropsSI('S','P',p_9pp,'Q',x_9pp,'Water')
    e_9pp = exergy(h_9pp,s_9pp)

    # State 9:
    p_9 = p_5
    T_9 = T_10p        # = T_10p = T_10pp = T_6
    h_9 = CP.PropsSI('H','P',p_9,'T',T_9,'Water')
    s_9 = CP.PropsSI('S','P',p_9,'T',T_9,'Water')
    x_9 = CP.PropsSI('Q','P',p_9,'T',T_9,'Water')
    e_9 = exergy(h_9,s_9)

    # State 6 -> LP inlet:
    p_6 = p_6  # = p_8
    T_6 = T_6  # = T_9 = T_10p = T_10pp
    h_6 = CP.PropsSI('H','P',p_6,'T',T_6,'Water')
    s_6 = CP.PropsSI('S','P',p_6,'T',T_6,'Water')
    x_6 = CP.PropsSI('Q','P',p_6,'T',T_6,'Water')
    e_6 = exergy(h_6,s_6)

    # State 8 -> Sup LP (6->8 isobar):
    p_8 = p_6
    T_8 = T_9p          # = T_9p = T_9pp
    h_8 = CP.PropsSI('H','P',p_8,'T',T_8,'Water')
    s_8 = CP.PropsSI('S','P',p_8,'T',T_8,'Water')
    x_8 = CP.PropsSI('Q','P',p_8,'T',T_8,'Water')
    e_8 = exergy(h_8,s_8)

    # State 8' -> Inlet Evap LP (2->8'->8"->8 isobar):
    p_8p = p_8
    x_8p = 0 # Saturated liquid
    T_8p = CP.PropsSI('T','P',p_8p,'Q',x_8p,'Water')
    h_8p = CP.PropsSI('H','P',p_8p,'Q',x_8p,'Water')
    s_8p = CP.PropsSI('S','P',p_8p,'Q',x_8p,'Water')
    e_8p = exergy(h_8p,s_8p)

    # State 8" -> Outlet Evap LP (2->8'->8"->8 isobar):
    p_8pp = p_8p
    x_8pp = 1 # Saturated vapor
    T_8pp = CP.PropsSI('T','P',p_8pp,'Q',x_8pp,'Water')
    h_8pp = CP.PropsSI('H','P',p_8pp,'Q',x_8pp,'Water')
    s_8pp = CP.PropsSI('S','P',p_8pp,'Q',x_8pp,'Water')
    e_8pp = exergy(h_8pp,s_8pp)

    # State 7 -> LP outlet:
    p_7 = p_7
    x_7 = x_7
    h_7 = CP.PropsSI('H','P',p_7,'Q',x_7,'Water')
    s_7 = CP.PropsSI('S','P',p_7,'Q',x_7,'Water')
    T_7 = CP.PropsSI('T','P',p_7,'Q',x_7,'Water')
    e_7 = exergy(h_7,s_7)

    # State 1 -> Inlet pump, Outlet Condenser (7->1 isothermal):
    T_1 = T_7
    x_1 = 0 # Saturated liquid
    p_1 = CP.PropsSI('P','T',T_1,'Q',x_1,'Water')
    h_1 = CP.PropsSI('H','T',T_1,'Q',x_1,'Water')
    s_1 = CP.PropsSI('S','T',T_1,'Q',x_1,'Water')
    e_1 = exergy(h_1,s_1)

    # State 2 -> Outlet pump (1->2 isothermal), Inlet Eco LP (2->8'->8"->8 isobar):
    p_2 = p_8pp
    s_2is = s_1
    h_2is = CP.PropsSI('H','P',p_2,'S',s_2is,'Water')
    h_2 = h_1 + (h_2is-h_1)/eta_pump
    T_2 = CP.PropsSI('T','P',p_2,'H',h_2,'Water')
    s_2 = CP.PropsSI('S','P',p_2,'H',h_2,'Water')
    x_2 = CP.PropsSI('Q','P',p_2,'H',h_2,'Water')
    e_2 = exergy(h_2,s_2)

    # State 4 -> HP outlet:
    p_4 = p_5 # Separator before the reheater
    s_4is = s_3
    h_4is = CP.PropsSI('H','P',p_4,'S',s_4is,'Water')
    h_4 = h_3 - eta_is_HP*(h_3-h_4is)
    s_4 = CP.PropsSI('S','P',p_4,'H',h_4,'Water')
    T_4 = CP.PropsSI('T','P',p_4,'H',h_4,'Water')
    x_4 = CP.PropsSI('Q','P',p_4,'H',h_4,'Water')
    e_4 = exergy(h_4,s_4)

    # State 5g -> Chimney:
    p_5g = p_4g

    T_g_HP = T_10p+T_pinch
    cp_g_HP = getMeanCp(p_4g,p_5g,T_4g,T_g_HP, R_f, comp, flue_conc_mass)
    h_g_HP = h_4g + cp_g_HP*(T_g_HP-T_4g)

    T_g_IP = T_9p+T_pinch
    cp_g_IP = getMeanCp(p_5g,p_5g,T_g_HP,T_g_IP, R_f, comp, flue_conc_mass)
    h_g_IP = h_g_HP + cp_g_IP*(T_g_IP-T_g_HP)

    T_g_LP = T_8p+T_pinch
    cp_g_LP = getMeanCp(p_5g,p_5g,T_g_IP,T_g_LP, R_f, comp, flue_conc_mass)
    h_g_LP = h_g_IP + cp_g_LP*(T_g_LP-T_g_IP)

    A = np.array([[(h_8-h_8p), (h_9p-h_8p), (h_9p-h_8p)], [(h_6-h_8), (h_9-h_9p), (h_10p-h_9p)], [0, (h_5-h_9), (h_3-h_10p+h_5-h_4)]])
    B = np.array([(h_g_IP-h_g_LP), (h_g_HP-h_g_IP), abs(h_4g-h_g_HP)])*dotm_g
    X = np.linalg.solve(A, B)

    dotm_vLP = X[0]
    dotm_vIP = X[1]
    dotm_vHP = X[2]

    dotm_v = dotm_vLP + dotm_vIP + dotm_vHP

    print("m_a = ", dotm_a)
    print("m_g = ", dotm_g)
    print("m_CH4 = ", dotm_f)
    print("m_vLP = ", dotm_vLP)
    print("m_vIP = ", dotm_vIP)
    print("m_vHP = ", dotm_vHP)
    print(50*'-')


    h_5g = h_4g - (dotm_vLP/dotm_g)*(h_6-h_2) - (dotm_vIP + dotm_vHP)*(h_5-h_2)/dotm_g
    T_5g = CP.PropsSI('T','P',p_5g,'H',h_5g*gas_prop[3],'Water')
    cp_5g = (h_5g-h_4g)/(T_5g-T_4g)
    s_5g = s_4g + cp_5g*np.log(T_5g/T_4g) - R_f*np.log(p_5g/p_4g)
    e_5g = (h_5g - h_1g) - T_1g*(s_5g - s_1g)

    e_exh = ((1 + 1/lamb_ma1)*e_5g - e_1g)*dotm_a


    ## Efficiencies calculations
    # ==========================
    # GAS TURBINE eta
    eta_cycleng = ((1+(1/lamb_ma1))*(h_3g - h_4g) - (h_2g - h_1g))/((1+(1/lamb_ma1))*h_3g - h_2g); #Energetic efficiency of the cycle
    eta_toteng = eta_mecg * eta_cycleng; #Total energetic efficiency = P_eg/(dotm_f*LHV)

    eta_cyclexg = ((1+(1/lamb_ma1))*(h_3g - h_4g) - (h_2g - h_1g))/((1+(1/lamb_ma1))*e_3g - e_2g); #Exergetic efficiency of the cycle
    eta_rotexg = ((1+(1/lamb_ma1))*(h_3g - h_4g) - (h_2g - h_1g))/((1+(1/lamb_ma1))*(e_3g - e_4g) - (e_2g - e_1g)); #Exergetic efficiency of the rotor assembly
    eta_combexg = (1/f)*((1+(1/lamb_ma1))*e_3g - e_2g)/((1+(1/lamb_ma1))*h_3g - h_2g); #Exergetic efficiency of the combustion
    eta_totexg = eta_mecg * eta_cyclexg * eta_combexg; #Total exergetic efficiency

    # STEAM TURBINE eta
    Wm_tot = (h_3-h_4)+(h_5-h_6)+(h_6-h_7)-(h_2-h_1)#-(h_9p-h_8p)-(h_10p-h_9p)
    eta_cyclenv = Wm_tot/(h_3+h_5+h_6-h_4-h_8p-h_9p-h_10p)
    eta_totenv = P_es/(dotm_g*(h_4g-h_5g))

    # GST eta
    eta_cyclen = eta_cycleng + eta_cyclenv - eta_cycleng*eta_cyclenv
    eps_exh = 0                                 # Losses at the chimney neglected (p.167 english book)
    eps_p = .01                                 # "we assume there are no unburnt residues etc" -p63
    eta_gen = 1-eps_p-eps_exh
    eta_toten = eta_toteng + eta_totenv*(1 - eta_toteng - eps_exh)


    ## Losses
    # =======

    emT = e_3-e_7
    WmT = h_3-h_7
    emTg = e_3g-e_4g
    WmTg = h_3g-h_4g

    emP = e_2-e_1
    WmP = h_2-h_1
    emPg = e_2g-e_1g
    WmPg = h_2g-h_1g

    loss_mec = k_mec*dotm_a*( (1 + 1/lamb_ma1)*(h_3g - h_4g) + (h_2g - h_1g) ) + (dotm_vLP + dotm_vIP + dotm_vHP)*Wm_tot - P_es #12.9e6
    loss_cond = dotm_v*(h_7 - h_1) #230.5e6
    loss_chimney = dotm_g*h_5g - dotm_a*h_1g #47.2e6

    loss_condex = dotm_v*(e_7 - e_1) #15.6e6
    loss_turbex = dotm_v*(emT - WmT) + dotm_g*(emTg - WmTg)
    loss_pumpex = dotm_v*(WmP - emP) + dotm_a*(WmPg - emPg)
    loss_rotex = loss_pumpex + loss_turbex #55.2e6
    loss_transex = dotm_g * (e_4g-e_5g) - dotm_vLP*( (e_6-e_2) ) - dotm_vIP*(e_5-e_2) - dotm_vHP*(e_3-e_2) #18.5e6
    loss_combex = dotm_a*(e_2g + e_c/lamb_ma1 - (1 + 1/lamb_ma1)*e_3g) #210.8e6
    loss_chemex = dotm_g*e_5g - dotm_a*e_1g #4.9e6
    loss_totex = loss_mec+loss_rotex+loss_combex+loss_chemex+loss_transex+loss_condex


    ## Figures
    # ========

    # 1st figure : Energetic balance
    fig_pie_en = plt.figure(1)
    labels = 'GT effective power \n'+ '%.1f'%(P_eg*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Condensor loss \n'+'%.1f'%(loss_cond*1e-6)+' MW', 'Chimney losses \n'+'%.1f'%(loss_chimney*1e-6)+' MW', 'ST effective power \n'+ '%.1f'%(P_es*1e-6)+' MW'
    sizes = [P_eg*1e-6, loss_mec*1e-6, loss_cond*1e-6, loss_chimney*1e-6, P_es*1e-6]
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=150)
    plt.axis('equal')
    plt.title("Primary power " + "%.1f" %(LHV*dotm_f*1e-6)+ " MW")


    # 2nd figure : Exergetic balance
    fig_pie_ex = plt.figure(2)
    labels = ['GT effective power \n'+ '%.1f'%(P_eg*1e-6)+' MW',    'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW',    'Condenser losses\n'+'%.1f'%(loss_condex*1e-6)+' MW',    'Rotor unit \n irreversibilities \n'+'%.1f'%(loss_rotex*1e-6)+' MW', 'Heat transfer irreversibilities \n'+'%.1f'%(loss_transex*1e-6)+' MW',    'Combustion \n irreversibilities \n'+'%.1f'%(loss_combex*1e-6)+' MW','Chimney losses \n'+'%.1f'%(loss_chemex*1e-6)+' MW', 'ST effective power \n'+ '%.1f'%(P_es*1e-6)+' MW']
    sizes = [P_eg*1e-6, loss_mec*1e-6, loss_condex*1e-6, loss_rotex*1e-6, loss_transex*1e-6,loss_combex*1e-6,loss_chemex*1e-6, P_es*1e-6]
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=65)
    plt.axis('equal')
    plt.title("Primary exergy flux " + "%.1f" %(e_c*dotm_f*1e-6) + " MW")
    
    # 3rd figure : T-s diagram
    fig_Ts_diagram = plt.figure(3)
    #fig_Ts_diagram = PropertyPlot('water', 'Ts', unit_system='EUR')
    #fig_Ts_diagram.calc_isolines(CoolProp.iQ, num=11)
    #fig_Ts_diagram.set_axis_limits([0., 9, 0, T_max+100-273.15])
    plt.grid(True)
    
    plt.title('T-s diagram of the cycle')
    plt.xlabel("s $[kJ/kg/K]$")
    plt.ylabel("t $[°C]$")
    
    
    # 4th figure: h-s diagram
    fig_hs_diagram = plt.figure(4)
    #fig_hs_diagram = PropertyPlot('water', 'HS', unit_system='EUR')
    #fig_hs_diagram.calc_isolines(CoolProp.iQ, num=11)
    #fig_hs_diagram.set_axis_limits([0., 9, 0, 4000])
    plt.grid(True)
    
    plt.title('h-s diagram of the cycle')
    plt.xlabel("s $[kJ/kg/K]$")
    plt.ylabel("h $[kJ/kg]$")
    
    # 5th figure: heat exchange diagram
    Q_4g = 0
    Q_g = (h_4g-h_5g)*dotm_g/dotm_vHP # /!\ Not right value!
    Q_3 = 0
    Q_Eco_LP = (h_8p-h_2)*dotm_v/dotm_vHP
    Q_Evap_LP = (h_8pp-h_8p)*dotm_vLP/dotm_vHP
    Q_Eco_IP = (h_9p-h_8p)*(dotm_vHP+dotm_vIP)/dotm_vHP
    Q_Sup_LP_LT = (h_6-h_8)*dotm_vLP/dotm_vHP
    Q_Evap_IP = (h_9pp-h_9p)*dotm_vIP/dotm_vHP
    Q_Sup_LP_HT = (h_8-h_8pp)*dotm_vLP/dotm_vHP
    Q_Sup_IP = (h_9-h_9pp)*dotm_vIP/dotm_vHP
    Q_Eco_HP = (h_10p-h_9p)*dotm_vHP/dotm_vHP
    Q_Evap_HP = (h_10pp-h_10p)*dotm_vHP/dotm_vHP
    Q_Sup_HP = (h_3-h_10pp)*dotm_vHP/dotm_vHP
    Q_Reh = (h_5-h_9)*dotm_vIP/dotm_vHP + (h_5-h_4)*dotm_vHP/dotm_vHP
    
    fig_heat_exchange = plt.figure(5)
    plt.grid(True)
    plt.plot([Q_4g, (Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP+Q_Sup_LP_LT+Q_Eco_IP+Q_Evap_LP+Q_Eco_LP)*1e-3], [T_4g-273.15, T_5g-273.15], 'r')
    plt.plot([Q_3, (Q_Sup_HP+Q_Reh)*1e-3], [T_3-273.15, T_10pp-273.15], 'b')
    plt.plot([(Q_Sup_HP+Q_Reh)*1e-3, (Q_Sup_HP+Q_Reh+Q_Evap_HP)*1e-3], [T_10pp-273.15, T_10p-273.15], 'b')
    plt.plot([(Q_Sup_HP+Q_Reh+Q_Evap_HP)*1e-3, (Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP)*1e-3], [T_10p-273.15, T_9pp-273.15], 'b')
    plt.plot([(Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP)*1e-3, (Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP)*1e-3], [T_9pp-273.15, T_9p-273.15], 'b')
    plt.plot([(Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP)*1e-3, (Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP+Q_Sup_LP_LT+Q_Eco_IP)*1e-3], [T_9p-273.15, T_8pp-273.15], 'b')
    plt.plot([(Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP+Q_Sup_LP_LT+Q_Eco_IP)*1e-3, (Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP+Q_Sup_LP_LT+Q_Eco_IP+Q_Evap_LP)*1e-3], [T_8pp-273.15, T_8p-273.15], 'b')
    plt.plot([(Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP+Q_Sup_LP_LT+Q_Eco_IP+Q_Evap_LP)*1e-3, (Q_Sup_HP+Q_Reh+Q_Evap_HP+Q_Sup_LP_HT+Q_Sup_IP+Q_Eco_HP+Q_Evap_IP+Q_Sup_LP_LT+Q_Eco_IP+Q_Evap_LP+Q_Eco_LP)*1e-3], [T_8p-273.15, T_2-273.15], 'b')
    
    plt.title('Recovery boiler heat exchange')
    plt.xlabel("Q $[kJ/kg_{vHP}]$")
    plt.ylabel("t $[°C]$")
    

    if display:
        plt.show()


    # Process output variables - do not modify---------------------------------
    p = (p_1g,p_2g,p_3g,p_4g,p_5g,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_8p,p_8pp,p_9,p_9p,p_9pp,p_10p,p_10pp)
    T = (T_1g,T_2g,T_3g,T_4g,T_5g,T_1,T_2,T_3,T_4,T_5,T_6,T_7,T_8,T_8p,T_8pp,T_9,T_9p,T_9pp,T_10p,T_10pp)
    s = (s_1g,s_2g,s_3g,s_4g,s_5g,s_1,s_2,s_3,s_4,s_5,s_6,s_7,s_8,s_8p,s_8pp,s_9,s_9p,s_9pp,s_10p,s_10pp)
    h = (h_1g,h_2g,h_3g,h_4g,h_5g,h_1,h_2,h_3,h_4,h_5,h_6,h_7,h_8,h_8p,h_8pp,h_9,h_9p,h_9pp,h_10p,h_10pp)
    e = (e_1g,e_2g,e_3g,e_4g,e_5g,e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_8p,e_8pp,e_9,e_9p,e_9pp,e_10p,e_10pp)
    x = (x_1g,x_2g,x_3g,x_4g,x_5g,x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_8p,x_8pp,x_9,x_9p,x_9pp,x_10p,x_10pp)
    DAT = (p,T,h,s,e,x)
    COMBUSTION = (LHV,e_c,excess_air,cp_gas,gas_prop)
    MASSFLOW = (dotm_a,dotm_f,dotm_g,dotm_v, dotm_vLP, dotm_vIP, dotm_vHP)
    DATEN = (loss_mec,loss_cond,loss_chimney)
    DATEX = (loss_mec,loss_rotex,loss_combex,loss_chemex,loss_transex,loss_totex,loss_condex)
    FIG = (fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram,fig_heat_exchange)
    out = (ETA,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG)
    return out
