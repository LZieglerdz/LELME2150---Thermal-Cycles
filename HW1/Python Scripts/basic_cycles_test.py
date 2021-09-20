"""
LELME2150 - Thermal cycles
Homework 1 - Basic cycles

Test code for your functions

Code execution can take a few seconds...

@author: Antoine Laterre
@date: September 10, 2021
"""

#
#===IMPORT PACKAGES============================================================
#

import numpy as np
import matplotlib.pyplot as plt
from basic_cycles_group_xx import gas_turbine
from basic_cycles_group_xx import steam_turbine

#
#===PLOT FUNCTIONS - IMPLEMENTED FOR YOU=======================================
#

def plot_eta_GT(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t):
    p_2_test, eta_p_2_test = np.linspace(2e+5,100e+5,100), np.zeros(100)    
    T_3_test, eta_T_3_test = np.linspace(873.15,1673.15,100), np.zeros(100)    
    eta_pi_test,  eta_eta_pi_test = np.linspace(0.7,1,100), np.zeros(100)    
    for i in range(100):
        p, T, s, h, eta_en = gas_turbine(p_1, T_1, p_2_test[i], p_2_test[i], T_3, eta_pi, eta_mec_c, eta_mec_t)
        eta_p_2_test[i] = eta_en
        p, T, s, h, eta_en = gas_turbine(p_1, T_1, p_2, p_3, T_3_test[i], eta_pi, eta_mec_c, eta_mec_t)
        eta_T_3_test[i] = eta_en
        p, T, s, h, eta_en = gas_turbine(p_1, T_1, p_2, p_3, T_3, eta_pi_test[i], eta_mec_c, eta_mec_t)
        eta_eta_pi_test[i] = eta_en        
    p, T, s, h, eta_en = gas_turbine(p_1, T_1, p_2, p_3, T_3, eta_pi, eta_mec_c, eta_mec_t)
    fig = plt.figure(figsize=[15,5])
    axs = fig.subplots(1,3)
    axs[0].plot(p_2_test*1e-5,eta_p_2_test,'-',color='tab:red',linewidth=2)
    axs[0].plot(p_2*1e-5,eta_en,'ok',linewidth=2)
    axs[0].set_xlabel('Compressor outlet p [bar]',fontsize=13)
    axs[0].set_xticks([0,20,40,60,80,100])
    axs[0].set_ylabel('GT energy efficiency [-]',fontsize=13)
    axs[0].set_yticks([0.0,0.1,0.2,0.3,0.4,0.5])    
    axs[1].plot(T_3_test-273.15,eta_T_3_test,'-',color='tab:red',linewidth=2)
    axs[1].plot(T_3-273.15,eta_en,'ok',linewidth=2)
    axs[1].set_xlabel('Combustor outlet T [°C]',fontsize=13)
    axs[1].set_xticks([600,800,1000,1200,1400])
    axs[1].set_yticks([0.0,0.1,0.2,0.3,0.4,0.5])
    axs[2].plot(eta_pi_test,eta_eta_pi_test,'-',color='tab:red',linewidth=2)
    axs[2].plot(eta_pi,eta_en,'ok',linewidth=2)
    axs[2].set_xlabel('Polytropic efficiency [-]',fontsize=13)
    axs[2].set_xticks([0.7,0.8,0.9,1.0])  
    axs[2].set_yticks([0.0,0.1,0.2,0.3,0.4,0.5])
    fig.tight_layout()
    plt.show()
    return fig

def plot_eta_ST(T_1,p_3,T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump):
    p_3_test, eta_p_3_test = np.linspace(60e+5,160e+5,100), np.zeros(100)    
    T_3_test, eta_T_3_test = np.linspace(623.15,1023.15,100), np.zeros(100)    
    eta_is_test,  eta_eta_is_test = np.linspace(0.7,1,100), np.zeros(100)    
    for i in range(100):
        p, T, s, h, x, eta_en, dot_m_f = steam_turbine(T_1,p_3_test[i],T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump)
        eta_p_3_test[i] = eta_en
        p, T, s, h, x, eta_en, dot_m_f = steam_turbine(T_1,p_3,T_3_test[i],eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump)
        eta_T_3_test[i] = eta_en
        p, T, s, h, x, eta_en, dot_m_f = steam_turbine(T_1,p_3,T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_test[i],eta_pump)
        eta_eta_is_test[i] = eta_en        
    p, T, s, h, x, eta_en, dot_m_f = steam_turbine(T_1,p_3,T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump)
    fig = plt.figure(figsize=[15,5])
    axs = fig.subplots(1,3)
    axs[0].plot(p_3_test*1e-5,eta_p_3_test,'-',color='tab:red',linewidth=2)
    axs[0].plot(p_3*1e-5,eta_en,'ok',linewidth=2)
    axs[0].set_xlabel('Boiler outlet p [bar]',fontsize=13)
    axs[0].set_xticks([60, 80, 100, 120, 140, 160])
    axs[0].set_ylabel('ST energy efficiency [-]',fontsize=13)
    axs[0].set_yticks([0.1,0.2,0.3])    
    axs[1].plot(T_3_test-273.15,eta_T_3_test,'-',color='tab:red',linewidth=2)
    axs[1].plot(T_3-273.15,eta_en,'ok',linewidth=2)
    axs[1].set_xlabel('Boiler outlet T [°C]',fontsize=13)
    axs[1].set_xticks([350,450,550,650,750])
    axs[1].set_yticks([0.1,0.2,0.3])
    axs[2].plot(eta_is_test,eta_eta_is_test,'-',color='tab:red',linewidth=2)
    axs[2].plot(eta_is_t,eta_en,'ok',linewidth=2)
    axs[2].set_xlabel('Isentropic efficiency [-]',fontsize=13)
    axs[2].set_xticks([0.7,0.8,0.9,1.0])  
    axs[2].set_yticks([0.1,0.2,0.3])
    fig.tight_layout()
    plt.show()
    return fig

#
#===BRAYTON CYCLE==============================================================
#

p_1, T_1 = 1e+5, 293.15 # [Pa], [K]
p_2 = 17.8e+5 # [Pa]
p_3, T_3 = 17.8e+5, 1273.15 # [Pa], [K]
eta_pi = 0.90 # [-]
eta_mec_c, eta_mec_t = 0.98, 0.98 # [-], [-]

p,T,s,h,eta_en = gas_turbine(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t)
fig_GT = plot_eta_GT(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t)

#
#===RANKINE CYCLE==============================================================
#

T_1 = 303.15 # [K]
p_3, T_3 = 100e+5, 813.15 # [Pa], [K]
eta_is_t = 0.85 # [-]
eta_mec_t = 0.985 # [-]
eta_pump = 0.8 # [-]
eta_gen = 0.6 # [-]
LHV = 25e+6 # [J/kg]
P_el = 95e+6 # [W]

p,T,s,h,x,eta_en,dot_m_f = steam_turbine(T_1,p_3,T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump)
fig_ST = plot_eta_ST(T_1,p_3,T_3,eta_gen,LHV,P_el,eta_mec_t,eta_is_t,eta_pump)