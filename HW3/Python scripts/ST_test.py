"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Test code steam turbine function.
Put this code in the same folder as the ST_group_xx.py file.

@author: Antoine Laterre
@date: October 22, 2021
"""

#
#===FILL YOUR GROUP N° HERE AND DO NOT MODIFY BELOW============================
#

GROUPS =  ['xx'] # Fill with your group number

#
#===IMPORT PACKAGES============================================================
#

import importlib

#
#===TEST CODE==================================================================
#

for gr in GROUPS:
    try:
        path = 'ST_group_' + gr;
        print('\n '+path)
        steam_turbine = importlib.import_module(path).steam_turbine
        P_e = 288e+6
        eta_turb = 0.92,0.90
        comb = 973.15,1.05,0.05,0.8
        options = 310e+5,70e+5,1e+5,288.15,838.15,305.15,353.15,0,0,0,421.85,0.88,comb,0.98,0.85,eta_turb
        display = True
        try:
            (ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG) = steam_turbine(P_e,options,display)
            (LHV,e_c,excess_air,cp_gas,gas_prop) = COMBUSTION
            (dotm_a,dotm_f,dotm_g,dotm_v) = MASSFLOW
            (eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex) = ETA
            (loss_mec,loss_gen,loss_cond) = DATEN
            (loss_mec,loss_rotex,loss_combex,loss_chemex,loss_transex,loss_totex,loss_condex) = DATEX
            (x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII) = XMASSFLOW
            (fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram) = FIG
            (p,T,h,s,e,x) = DAT
            (p_1,p_2,p_3,p_4,p_5,
              p_6,p_6I,p_6II,p_6III,p_6IV,p_6V,p_6VI,p_6VII,p_6VIII,
              p_7,p_7I,p_7II,p_7III,p_7IV,p_7V,p_7VI,p_7VII,p_7VIII,
              p_8,
              p_9,p_9I,p_9II,p_9III,p_9IV,p_9V,p_9VI,p_9VII,p_9VIII) = p
            (T_1,T_2,T_3,T_4,T_5,
              T_6,T_6I,T_6II,T_6III,T_6IV,T_6V,T_6VI,T_6VII,T_6VIII,
              T_7,T_7I,T_7II,T_7III,T_7IV,T_7V,T_7VI,T_7VII,T_7VIII,
              T_8,
              T_9,T_9I,T_9II,T_9III,T_9IV,T_9V,T_9VI,T_9VII,T_9VIII) = T
            (s_1,s_2,s_3,s_4,s_5,
              s_6,s_6I,s_6II,s_6III,s_6IV,s_6V,s_6VI,s_6VII,s_6VIII,
              s_7,s_7I,s_7II,s_7III,s_7IV,s_7V,s_7VI,s_7VII,s_7VIII,
              s_8,
              s_9,s_9I,s_9II,s_9III,s_9IV,s_9V,s_9VI,s_9VII,s_9VIII) = s
            (h_1,h_2,h_3,h_4,h_5,
              h_6,h_6I,h_6II,h_6III,h_6IV,h_6V,h_6VI,h_6VII,h_6VIII,
              h_7,h_7I,h_7II,h_7III,h_7IV,h_7V,h_7VI,h_7VII,h_7VIII,
              h_8,
              h_9,h_9I,h_9II,h_9III,h_9IV,h_9V,h_9VI,h_9VII,h_9VIII) = h
            (e_1,e_2,e_3,e_4,e_5,
              e_6,e_6I,e_6II,e_6III,e_6IV,e_6V,e_6VI,e_6VII,e_6VIII,
              e_7,e_7I,e_7II,e_7III,e_7IV,e_7V,e_7VI,e_7VII,e_7VIII,
              e_8,
              e_9,e_9I,e_9II,e_9III,e_9IV,e_9V,e_9VI,e_9VII,e_9VIII) = e
            (x_1,x_2,x_3,x_4,x_5,
              x_6,x_6I,x_6II,x_6III,x_6IV,x_6V,x_6VI,x_6VII,x_6VIII,
              x_7,x_7I,x_7II,x_7III,x_7IV,x_7V,x_7VI,x_7VII,x_7VIII,
              x_8,
              x_9,x_9I,x_9II,x_9III,x_9IV,x_9V,x_9VI,x_9VII,x_9VIII) = x        
            print('\n Basic ST passed: your code will be graded!')
        except:
            print('\n Basic ST failed: please correct your code!')
    except:
        print('\n Basic ST failed: please correct your code!')