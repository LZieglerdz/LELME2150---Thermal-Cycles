"""
LELME2150 - Thermal cycles
Homework 2 - Exergy and combustion

Signature of the final gas turbine function

@author: Antoine Laterre
@date: September 22, 2021
"""

#
#===IMPORT PACKAGES============================================================
#

import CoolProp.CoolProp as CP
import numpy as np
from thermochem import janaf
db = janaf.Janafdb();
from scipy.integrate import quad;
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
# import re      #try to use regex to get composition of fuel from chemical formula: C4H10 -> C:4, H:10;

#
#===GLOBAL VARIABLES===========================================================
#
T_S = 273.153         # [K]
p_S = 1e5             # [Pa]
fuel = 'CH4'
# LHV = 50150e3         # J/kg_CH4

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
Mm_fuel = CP.PropsSI("MOLARMASS", fuel)        # kg/mol
Mm_C = Mm_CO2-Mm_O2        # kg/mol
Mm_H = .5*Mm_H2         # kg/mol
Mm_O = .5*Mm_O2      # kg/mol

comp = ["O2","N2","CO2","H2O"]
air_conc = [O2_conc,N2_conc,CO2_conc,H2O_conc]
Mm   = [Mm_O2,Mm_N2,Mm_CO2,Mm_H2O]
Mm_air = np.dot(Mm, air_conc)

def get_ma1(x, y): #Stoechiometric air-to-fuel ratio [-]
    return (Mm_O2+3.76*Mm_N2)*(1+(y-2*x)*0.25)/(Mm_C+Mm_H*y+Mm_O*x);

def get_LHV(x, y): # Lower Heating Value for solid fuels CHyOx
    LHV_mol = 393400 + 102250*y - x*(111000 + 102250*y)/(1 + y/2) # [kJ/kmol]
    LHV = LHV_mol/(12 + y + 16*x) # [kJ/kg]
    return LHV

def CombEx(x, y): # CHyOx # e_c and Cp of different fuels (from slides of LMECA2150)
    if x==0: # Combustible of type CHy
        if y==0: # C
            e_c = 34160 # [kJ/kg]
            Cp = 0.86667
        if y==1.8: # CH1.8
            e_c = 45710 # [kJ/kg]
            Cp = 1.10434
        if y==4: # CH4
            e_c = 52215 # [kJ/kg]
            Cp = 2.2005
    elif y==0: # Combustible of type COx
        if x==1: # CO
            e_c = 9845 # [kJ/kg]
            Cp = 1.71176
    else:
        m_water = (y/2)*18/1000; # [kg]
        HHV = get_LHV(x,y) + 2375*m_water; # [kJ/kg]
        e_c = HHV - 5350; # [kJ/kg] T0*S = 5350 (S is the carbon entropy at standard conditions)
        Cp = 0.86667;
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

def get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3,iter, lam_est):
    #  Combustion supposée complète
    #  C_zH_yO_x + w(O_2 + 3.76 N_2) --> a_0 O_2 + a_2 CO_2 + b_2 H_2O + 3.76w N_2
    # comb stoechiometric -> w = z + (y-2*x)/4
    # a_2 = z
    # 2*b_2 = y
    # 2*a_0 + 2*a_2 + b_2 = x + 2*w
    w = lam_est*(z + (y-2*x)/4)
    a_0 = (x + 2*w - (y/2) - 2*z)/2
    Mm_f, flue_conc_mol, R_f = FlueGasFrac(x, y, lam_est)
    flue_conc_mass = np.multiply(flue_conc_mol, [Mm_O2, Mm_N2, Mm_CO2, Mm_H2O])/Mm_f

    cp_3S = getMeanCp(p_S, p_3, T_S, T_3, R_f, comp, flue_conc_mass)
    cp_32 = getMeanCp(p_2, p_3, T_2, T_3, R_f, comp, flue_conc_mass)

    lam = (get_LHV(x,y)*1e3-cp_3S*(T_3-T_S))/ (get_ma1(x,y)*cp_32*(T_3-T_2))

    if iter == 0:
        print('The function get_lambda did not converge.')
        return(-1)
    if np.abs(lam-lam_est) <= 1e-6 :
        return(lam)
    else:
        return( get_lambda(fuel, z, y, x, T_2, T_3, p_2, p_3, iter-1, lam) )

#
#===MEAN CP====================================================================
#

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
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

# IN:
#   - P_e [W]: net power output
#   - options: tuple containing the parametric values of the GT
#       o p_1 [Pa], T_1 [K]: inlet pressure and temperature (ambient)
#       o r [-], eta_pi_c [-]: compression ratio, compressor polytropic efficiency
#       o T_3 [K], k_cc [-]: combustor outlet temperature and pressure losses coeffcient
#       o eta_pi_t [-]: turbine polytropic efficiency
#       o k_mec [-]: shaft losses
#   - display: bool to choose to plot the T-s & h-s diagrams and the energy and exergy pie charts (True or False)
# OUT: tuple containing
#   - DAT: tuple containing the cycle state data
#      o p [Pa]: tuple containing the pressure at each state
#      o T [T]: tuple containing the temperature at each state
#      o h [J/kg]: tuple containing the enthalpy at each state
#      o s [J/kg/K]: tuple containing the entropy at each state
#      o e [J/kg]: tuple containing the exergy at each state
#   - COMBUSTION: tuple containing the combustion parameters
#      o LHV [J/kg]: the fuel Lower Heating Value
#      o e_c [J/kg]: the fuel exergy
#      o excess_air [-]: the excess air of the combustion
#      o cp_gas [J/kg/K]: the flue gas specific heat capacity at the combustor outlet
#      o gas_prop: list containing the proportion of ['CO2','H2O','N2','O2'] in the flue gas respectively
#   - MASSFLOW: tuple containing the massflow rates
#      o dot_m_a [kg/s]: mass flow rate of air
#      o dot_m_f [kg/s]: mass flow rate of fuel
#      o dot_m_g [kg/s]: mass flow rate of flue gas
#   - ETA: tuple containing the efficiencies (see text book pp. 113-125)
#       o eta_cyclen [-]: cycle energy efficiency
#       o eta_toten [-]: overall energy efficiency
#       o eta_cyclex [-]: cycle exergy efficiency
#       o eta_totex [-]: overall exergy efficiency
#       o eta_rotex [-]: compressor-turbine exergy efficiency
#       o eta_combex [-]: combustion exergy efficiency
#   - DATEN: tuple containing the energy losses
#       o loss_mec [W]: mechanical energy losses
#       o loss_ech [W]: exhaust energy losses
#   - DATEX: tuple containing the exergy losses
#       o loss_mec [W]: mechanical energy losses
#       o loss_rotex [W]: compressor-turbine exergy losses
#       o loss_combex [W]: combustion exergy losses
#       o loss_echex [W]: exhaust exergy losses
#   - FIG: tuple containing the figures to be diplayed ((0,0,0,0) if display = False)
#       o fig_pie_en: pie chart of energy losses
#       o fig_pie_ex: pie chart of exergy losses
#       o fig_Ts_diagram: T-s diagram of the GT cycle
#       o fig_hs_diagram: h-s diagram of the GT cycle
# REMARK:
#   - You can adapt the numerical application of the text book
#   in page pp. 125 to check your results.
# - As improvements, other options can be added (not expected before HW4):
#       o air [n.a.], air_prop [-]: air composition (types: list of str, list of float)
#       o NTU [-]: number of transfer unit in case of a recuperator
#       o fuel [n.a.]: type of fuel  (type: str), default is 'methane'
def gas_turbine(P_e,options,display):
    # Process input variables--------------------------------------------------
    p_1,T_1,r,eta_pi_c,T_3,k_cc,eta_pi_t,k_mec = options

    # Replace with your model--------------------------------------------------
    p_1, p_2, p_3, p_4 = p_1, 0, 0, 0
    T_1, T_2, T_3, T_4 = T_1, 0, T_3, 0
    h_1, h_2, h_3, h_4 = 0, 0, 0, 0
    s_1, s_2, s_3, s_4 = 0, 0, 0, 0
    e_1, e_2, e_3, e_4 = 0, 0, 0, 0
    LHV,e_c,excess_air,cp_gas,gas_prop = 0, 0, 0, 0, 0
    dotm_a,dotm_f,dotm_g = 0, 0, 0
    eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex = 0, 0, 0, 0, 0, 0
    loss_mec,loss_ech = 0, 0
    loss_mec,loss_rotex,loss_combex,loss_echex = 0, 0, 0, 0
    fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram = 0, 0, 0, 0

    LHV = get_LHV(0,4); #Lower heating value of CH4 [kJ/kg_CH4]
    e_c = CombEx(0, 4)[1]
    f = e_c/LHV; #e_comb/LHV
    ma1 = get_ma1(0,4);
    Cp_c = CombEx(0, 4)[0]

    R_Star = 0
    for i in range(len(comp)):
        R_Star += air_conc[i]*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, comp[i])/Mm[i] # [J/kgK]

    #set references state
    for i in comp:
        Dmolar = CP.PropsSI("Dmolar", "T", T_S, "P", p_1, i)
        CP.set_reference_state(i, T_S, Dmolar, 0, 0)


    # State 1 -- 4->1: Isobar Heat Rejection
    s_1 = .79*CP.PropsSI('S','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('S','T', T_1,'P', p_1, "O2")
    h_1 = .79*CP.PropsSI('H','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('H','T', T_1,'P', p_1, "O2")
    e_1 = (h_1 - h_1) - T_1*(s_1 - s_1)

    # State 2 -- 1->2: Polytropic Compression
    p_2 = r*p_1
    T_2 = getPolytropicTemp(p_1, p_2, T_1, T_1, R_Star, eta_pi_c, 100, comp, air_conc)
    h_2 = h_1 + getMeanCp(p_1,p_2,T_1,T_2, R_Star,comp, air_conc)*(T_2-T_1)
    s_2 = s_1 + getMeanCp(p_1,p_2,T_1,T_2, R_Star, comp, air_conc)*np.log(T_2/T_1) - R_Star*np.log(p_2/p_1)
    e_2 = (h_2 - h_1) - T_1*(s_2 - s_1)

    # State 3 -- 2->3: Isobar Combustion
    p_3 = k_cc*p_2
    lamb = get_lambda('CH4', 1, 4, 0, T_2, T_3, p_2, p_3, 200, 1)
    print("lambda: %.2f" %lamb)

    Mm_f, flue_conc_mol, R_f = FlueGasFrac(0, 4, lamb)
    flue_conc_mass = np.multiply(flue_conc_mol, [Mm_O2, Mm_N2, Mm_CO2, Mm_H2O])/Mm_f

    h_3 = h_2
    s_3 = s_2
    h_3 += getMeanCp(p_2,p_3,T_2,T_3, R_f, comp, flue_conc_mass)*(T_3-T_2)
    s_3 += getMeanCp(p_2,p_3,T_2,T_3, R_f, comp, flue_conc_mass)*np.log(T_3/T_2) - R_f*np.log(p_3/p_2)
    e_3 = (h_3 - h_1) - T_1*(s_3 - s_1)

    lamb_ma1 = lamb*ma1;

    # State 4 -- 3->4: Polytropic Expansion
    p_4 = p_1
    T_4 = getPolytropicTemp(p_3, p_4, T_3, T_3, R_f , 1/eta_pi_t, 100, comp, flue_conc_mass)

    h_4 = h_3 + getMeanCp(p_3,p_4,T_4,T_3, R_f, comp, flue_conc_mass)*(T_4-T_3)
    s_4 = s_3 + getMeanCp(p_3,p_4,T_3,T_4, R_f, comp, flue_conc_mass)*np.log(T_4/T_3) - R_f*np.log(p_4/p_3)
    e_4 = (h_4 - h_1) - T_1*(s_4 - s_1)

    ## Efficiencies calculations
    # ==========================

    eta_cyclen = ((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1))/((1+(1/lamb_ma1))*h_3 - h_2); #Energetic efficiency of the cycle
    eta_mec = 1 - k_mec*((1+(1/lamb_ma1))*(h_3 - h_4) + (h_2 - h_1))/((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1)); #Mechanical efficiency
    eta_toten = eta_mec * eta_cyclen; #Total energetic efficiency

    eta_cyclex = ((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1))/((1+(1/lamb_ma1))*e_3 - e_2); #Exergetic efficiency of the cycle
    eta_rotex = ((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1))/((1+(1/lamb_ma1))*(e_3 - e_4) - (e_2 - e_1)); #Exergetic efficiency of the rotor assembly
    eta_combex = (1/f)*((1+(1/lamb_ma1))*e_3 - e_2)/((1+(1/lamb_ma1))*h_3 - h_2); #Exergetic efficiency of the combustion
    eta_totex = eta_mec * eta_cyclex * eta_combex; #Total exergetic efficiency

    ## Energetic losses
    # =================

    flow_air = P_e/(eta_mec*((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1)));  #Air mass flow rate [kg_air/s]

    loss_mec = k_mec*((1+(1/lamb_ma1))*(h_3 - h_4) + (h_2 - h_1))*flow_air; #Mechanical losses [W]
    loss_ech = ((1+(1/lamb_ma1))*h_4 - h_1)*flow_air;                       #Exhaust gas losses [W]

    ## Mass flows
    # ===========

    dotm_a = flow_air; #Air mass flow rate [kg_air/s]
    dotm_f = flow_air/lamb_ma1; #Combustible mass flow rate [kg_comb/s]
    dotm_g = (1+(1/lamb_ma1))*flow_air;#Fluegas mass flow rate [kg_fluegas/s]

    excess_air = lamb;
    x_O2f, x_N2f, x_CO2f, x_H2Of = flue_conc_mol
    gas_prop = [x_CO2f, x_H2Of, x_N2f, x_O2f]
    cp_gas = getMeanCp(p_S, p_3, T_S, T_3, R_f, comp, flue_conc_mass) #cp_fluegas(T0,400); #Cp of fluegas at 400[K] [kJ/kg_fluegas]

    # gas_prop[0] = x_O2f*(1+(1/lamb_ma1))*flow_air;  #Mass flow rate of O2 in exhaust gas [kg_O2/s]
    # gas_prop[1] = x_N2f*(1+(1/lamb_ma1))*flow_air;  #Mass flow rate of N2 in exhaust gas [kg_N2/s]
    # gas_prop[2] = x_CO2f*(1+(1/lamb_ma1))*flow_air; #Mass flow rate of CO2 in exhaust gas [kg_CO2/s]
    # gas_prop[3] = x_H2Of*(1+(1/lamb_ma1))*flow_air; #Mass flow rate of H2O in exhaust gas [kg_H2O/s]


    ## Exergetic losses
    # =================

    # pb avec loss_combex
    perte_combex = dotm_g*(e_c - s_2 + s_3 + s_4 ) - (dotm_f*e_3-dotm_a*e_2)     #Ju
    # loss_combex = flow_air*(e_2 + f*LHV/lamb_ma1 - (1+(1/lamb_ma1))*e_3); #Combustion losses [W]
    print(e_2)
    print( f, LHV, lamb_ma1)
    print((1+(1/lamb_ma1))*e_3)
    loss_turbine = flow_air*(1+(1/lamb_ma1))*((e_3 - e_4) - (h_3 - h_4)); #Shaft losses at the turbine [W]
    loss_compressor = flow_air*((h_2 - h_1) - (e_2 - e_1));               #Shaft losses at the compressor [W]
    loss_rotex = (loss_compressor + loss_turbine);                        #Total shaft losses [W]
    loss_echex = flow_air*((1+(1/lamb_ma1))*e_4 - e_1);                   #Exhaust gas losses [W]

    ## Generate graph to export:
    # ==========================

    # My 1st figure: Pie chart energy
    fig_pie_en = plt.figure(1)
    labels = 'Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Exhaust losses MW \n'+'%.1f'%(loss_ech*1e-6)+' MW'
    sizes = [P_e*1e-6, loss_mec*1e-6, loss_ech*1e-6]
    colors = ['gold', 'yellowgreen', 'lightcoral']
    #explode = (0, 0, 0, 0)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')
    plt.title("Primary energy flux " + "%.1f" %(LHV*dotm_f*1e-6) + " MW")

    # # My 2nd figure: Pie chart exergy
    fig_pie_ex = plt.figure(2)
    labels = 'Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Exhaust losses \n'+'%.1f'%(loss_echex*1e-6)+' MW', 'Turbine & compressor \n irreversibilities \n'+'%.1f'%(loss_rotex*1e-6)+' MW', 'Combustion \n irreversibilities \n'+'%.1f'%(loss_combex*1e-6)+' MW'
    sizes = [P_e*1e-6, loss_mec*1e-6, loss_echex*1e-6, loss_rotex*1e-6, loss_combex*1e-6]
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'red']
    #explode = (0, 0, 0, 0)
    # plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=340)
    # plt.axis('equal')
    # plt.title("Primary exergy flux " + "%.1f" %(e_c*dotm_f*1e-6) + " MW")
    print(sizes)
    plt.show()

    # Process output variables - do not modify---------------------------------
    p = (p_1, p_2, p_3, p_4)
    T = (T_1, T_2, T_3, T_4)
    s = (s_1, s_2, s_3, s_4)
    h = (h_1, h_2, h_3, h_4)
    e = (e_1, e_2, e_3, e_4)
    DAT = (p,T,h,s,e)
    COMBUSTION = (LHV,e_c,excess_air,cp_gas,gas_prop)
    MASSFLOW = (dotm_a,dotm_f,dotm_g)
    ETA = (eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex)
    DATEN = (loss_mec,loss_ech)
    DATEX = (loss_mec,loss_rotex,loss_combex,loss_echex)
    FIG = (fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram)
    out = (ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG)
    return out
