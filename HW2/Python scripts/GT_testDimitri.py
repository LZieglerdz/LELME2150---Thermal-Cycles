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

#import warnings
#warnings.filterwarnings('ignore')


#
#=== MEAN CP ==========================================
#
def getMeanCp(p_in, p_out, T_in, T_out, R_Star):
    T_min = np.min([T_in,T_out])
    T_max = np.max([T_in,T_out])
    p_min = np.min([p_in,p_out])
    p_max = np.max([p_in,p_out])
    cp = .79*CP.PropsSI('CPMASS','T', T_min,'P', p_min, "N2") + .21*CP.PropsSI('CPMASS','T', T_min,'P', p_min, "O2")
    if T_min == T_max:
        return(cp)
    n=100
    rangeT = np.linspace(T_min,T_max,n)
    rangeP = np.ones(n)*p_min
    p = p_min
    for i in np.arange(1,n):
        gamma = cp/(cp - i*R_Star)
        p = p_min*(rangeT[i]/T_min)**((gamma-1)/gamma)
        cp += .79*CP.PropsSI('CPMASS','T', rangeT[i],'P', rangeP[i], "N2") + .21*CP.PropsSI('CPMASS','T', rangeT[i],'P', p, "O2")
    return(cp/n)

def getPolytropicTemp(p_in, p_out, T_in, T_out, R_Star, eta_pi, iter):
    if iter < 0:
        print("Function does not converge.")
        return(-1)
    T = T_in*(p_out/p_in)**(R_Star/(eta_pi*getMeanCp(p_in,p_out,T_in,T_out, R_Star)))
    if np.abs(T_out-T) < 1e-12:
        return(T)
    else:
        return(getPolytropicTemp(p_in, p_out, T_in, T, R_Star, eta_pi, iter-1))

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
    
    R_Star = .79*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, "N2")/CP.PropsSI('MOLAR_MASS','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('GAS_CONSTANT','T', T_1,'P', p_1, "O2")/CP.PropsSI('MOLAR_MASS','T', T_1,'P', p_1, "O2") # [J/kgK]

    #set references state
    DmolarN2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "N2")
    CP.set_reference_state("N2", T_1, DmolarN2, 0, 0)
    DmolarO2 = CP.PropsSI("Dmolar", "T", T_1, "P", p_1, "O2")
    CP.set_reference_state("O2", T_1, DmolarO2, 0, 0)

    # State 1 -- 4->1: Isobar Heat Rejection
    s_1 = .79*CP.PropsSI('S','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('S','T', T_1,'P', p_1, "O2")
    h_1 = .79*CP.PropsSI('H','T', T_1,'P', p_1, "N2") + .21*CP.PropsSI('H','T', T_1,'P', p_1, "O2")
    e_1 = 0; # Reference state

    # State 2 -- 1->2: Polytropic Compression
    p_2 = r*p_1;
    T_2 = getPolytropicTemp(p_1, p_2, T_1, T_1, R_Star, eta_pi_c, 100)
    h_2 = getMeanCp(p_1,p_2,T_1,T_2, R_Star)*(T_2-T_1)
    s_2 = s_1 + getMeanCp(p_1,p_2,T_1,T_2, R_Star)*np.log(T_2/T_1) - R_Star*np.log(p_2/p_1)
    e_2 = (h_2 - h_1) - T_1*(s_2 - s_1)

    # State 3 -- 2->3: Isobar Combustion
    p_3 = k_cc*p_2;
    
    Mm_O2 = CP.PropsSI('MOLARMASS', 'P', p_3, 'T', T_3, 'O2') # [kg/kmol]
    Mm_N2 = CP.PropsSI('MOLARMASS', 'P', p_3, 'T', T_3, 'N2')
    Mm_CO2 = CP.PropsSI('MOLARMASS', 'P', p_3, 'T', T_3, 'CO2')
    Mm_H2O = CP.PropsSI('MOLARMASS', 'P', p_3, 'T', T_3, 'H2O')
    Mm_CO = CP.PropsSI('MOLARMASS', 'P', p_3, 'T', T_3, 'CO')


    def ma1(x, y): #Stoechiometric air-to-fuel ratio [-]
        return (Mm_O2+3.76*Mm_N2)*(1+(y-2*x)*0.25)/(12.01+1.008*y+16*x);
    
    def LHV(x, y): # Lower Heating Value for solid fuels CHyOx
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
            HHV = LHV(x,y) + 2375*m_water; # [kJ/kg]
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
        R_f = janaf.R*1000/Mm_f; # Ideal gas constant (R*) for exhaust gas [J/kg/K]
       
        return [Mm_f, [x_O2f, x_N2f, x_CO2f, x_H2Of], R_f]
    
    LHV = LHV(0,4); #Lower heating value of CH4 [kJ/kg_CH4]
    e_c = CombEx(0, 4)[1]
    f = e_c/LHV; #e_comb/LHV
    ma1 = ma1(0,4);
    Cp_c = CombEx(0, 4)[0]
    
    """ Trouver :
        lamb = excess air coefficient
        h_fluegas = h_3
        s_fluegas = s_3
        e_fluegas = e_3
        T_fluegas = T_3 (maximal T ? 1600K ?)
        Cp_fluegas
    """


    lamb_ma1 = lamb*ma1;
    
    excess_air = lamb;
    cp_gas = 1 #cp_fluegas(T0,400); #Cp of fluegas at 400[K] [kJ/kg_fluegas]
    gas_prop[0] = m_O2f*(1+(1/lamb_ma1))*flow_air;  #Mass flow rate of O2 in exhaust gas [kg_O2/s]
    gas_prop[1] = m_N2f*(1+(1/lamb_ma1))*flow_air;  #Mass flow rate of N2 in exhaust gas [kg_N2/s]
    gas_prop[2] = m_CO2f*(1+(1/lamb_ma1))*flow_air; #Mass flow rate of CO2 in exhaust gas [kg_CO2/s]
    gas_prop[3] = m_H2Of*(1+(1/lamb_ma1))*flow_air; #Mass flow rate of H2O in exhaust gas [kg_H2O/s]
    
    
    s_3 = .79*CP.PropsSI('S','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('S','T', T_3,'P', p_3, "O2")
    h_3 = .79*CP.PropsSI('H','T', T_3,'P', p_3, "N2") + .21*CP.PropsSI('H','T', T_3,'P', p_3, "O2")

    # State 4 -- 3->4: Polytropic Expansion
    p_4 = p_1
    T_4 = getPolytropicTemp(p_3, p_4, T_3, T_3, R_Star, 1/eta_pi_t, 100)
    h_4 = h_3 + getMeanCp(p_3,p_4,T_4,T_3, R_Star)*(T_4-T_3)
    s_4 = s_3 + getMeanCp(p_3,p_4,T_3,T_4, R_Star)*np.log(T_4/T_3) - R_Star*np.log(p_4/p_3)

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
    
    ## Exergetic losses
    # =================
    
    loss_combex = flow_air*(e_2 + f*LHV/lamb_ma1 - (1+(1/lamb_ma1))*e_3); #Combustion losses [W]
    loss_turbine = flow_air*(1+(1/lamb_ma1))*((e_3 - e_4) - (h_3 - h_4)); #Shaft losses at the turbine [W]
    loss_compressor = flow_air*((h_2 - h_1) - (e_2 - e_1));               #Shaft losses at the compressor [W]
    loss_rotex = (loss_compressor + loss_turbine);                        #Total shaft losses [W]
    loss_echex = flow_air*((1+(1/lamb_ma1))*e_4 - e_1);                   #Exhaust gas losses [W]
    
    ## Mass flows
    # ===========
    
    dotm_a = flow_air; #Air mass flow rate [kg_air/s]
    dotm_f = flow_air/lamb_ma1; #Combustible mass flow rate [kg_comb/s]
    dotm_g = (1+(1/lamb_ma1))*flow_air;#Fluegas mass flow rate [kg_fluegas/s]
    
    ## Generate graph to export:
    # ==========================
   
    # My 1st figure: Pie chart energy
    fig_pie_en = plt.figure(1)
    labels = 'Effective Power \n'+ '%.1f'%(Pe*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Exhaust losses MW \n'+'%.1f'%(loss_ech*1e-6)+' MW'
    sizes = [Pe*1e-6, loss_mec*1e-6, loss_ech*1e-6]
    colors = ['gold', 'yellowgreen', 'lightcoral']
    #explode = (0, 0, 0, 0)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')
    plt.title("Primary energy flux " + "%.1f" %(LHV*dotm_f*1e-6) + " MW")
   
    # My 2nd figure: Pie chart exergy
    fig_pie_ex = plt.figure(2)
    labels = 'Effective Power \n'+ '%.1f'%(Pe*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Exhaust losses \n'+'%.1f'%(loss_echex*1e-6)+' MW', 'Turbine & compressor \n irreversibilities \n'+'%.1f'%(loss_rotex*1e-6)+' MW', 'Combustion \n irreversibilities \n'+'%.1f'%(loss_combex*1e-6)+' MW'
    sizes = [Pe*1e-6, loss_mec*1e-6, loss_echex*1e-6, loss_rotex*1e-6, loss_combex*1e-6]
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'red']
    #explode = (0, 0, 0, 0)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=340)
    plt.axis('equal')
    plt.title("Primary exergy flux " + "%.1f" %(e_c*dotm_f*1e-6) + " MW")
    
    
    
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

