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

GROUPS =  ['42'] # Fill with your group number

#
#===IMPORT PACKAGES============================================================
#

import importlib
from tabulate import tabulate

#
#===TEST CODE==================================================================
#
# gr = '42'
# path = 'ST_group_' + gr;
# print('\n '+path)
# steam_turbine = importlib.import_module(path).steam_turbine
# P_e = 288e+6
# eta_turb = 0.92,0.90
# comb = 973.15,1.05,0.05,0.8
# p_3 = 310e+5
# p_4 = 70e+5
# p_ref = 1e+5
# T_ref = 288.15
# T_max = 838.15
# T_cond_out = 305.15
# T_exhaust = 353.15
# T_pinch_sub = 0
# T_pinch_ex = 0
# T_pinch_cond = 1
# T_drum = 421.85
# x_6 = 0.88
# eta_mec = 0.98
# eta_pump = 0.85
# options = p_3,p_4,p_ref,T_ref,T_max,T_cond_out,T_exhaust,T_pinch_sub,T_pinch_ex,T_pinch_cond,T_drum,x_6,comb,eta_mec,eta_pump,eta_turb
# display = False
# steam_turbine(P_e,options,display)


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

    header  = ['States', 't\n°C', 'p\nkPa', 'x\n-', 'h\nkJ/kg', 's\nkJ/kgK', 'e\nkJ/kg']
    data    = [ ['1    ','%.2f'%(T[0]-273.15), '%.1f'%(p[0]*1e-3), '%.2f'%(x[0]), '%.1f'%(h[0]*1e-3), '%.3f'%(s[0]*1e-3), '%.1f'%(e[0]*1e-3)],
                ['2    ','%.2f'%(T[1]-273.15), '%.1f'%(p[1]*1e-3), '%.2f'%(x[1]), '%.1f'%(h[1]*1e-3), '%.3f'%(s[1]*1e-3), '%.1f'%(e[1]*1e-3)],
                ['3    ','%.2f'%(T[2]-273.15), '%.1f'%(p[2]*1e-3), '%.2f'%(x[2]), '%.1f'%(h[2]*1e-3), '%.3f'%(s[2]*1e-3), '%.1f'%(e[2]*1e-3)],
                ['4    ','%.2f'%(T[3]-273.15), '%.1f'%(p[3]*1e-3), '%.2f'%(x[3]), '%.1f'%(h[3]*1e-3), '%.3f'%(s[3]*1e-3), '%.1f'%(e[3]*1e-3)],
                ['6VIII','%.2f'%(T[13]-273.15), '%.1f'%(p[13]*1e-3), '%.2f'%(x[13]), '%.1f'%(h[13]*1e-3), '%.3f'%(s[13]*1e-3), '%.1f'%(e[13]*1e-3)],
                ['5    ','%.2f'%(T[4]-273.15), '%.1f'%(p[4]*1e-3), '%.2f'%(x[4]), '%.1f'%(h[4]*1e-3), '%.3f'%(s[4]*1e-3), '%.1f'%(e[4]*1e-3)],
                ['6VII ','%.2f'%(T[12]-273.15), '%.1f'%(p[12]*1e-3), '%.2f'%(x[12]), '%.1f'%(h[12]*1e-3), '%.3f'%(s[12]*1e-3), '%.1f'%(e[12]*1e-3)],
                ['6VI  ','%.2f'%(T[11]-273.15), '%.1f'%(p[11]*1e-3), '%.2f'%(x[11]), '%.1f'%(h[11]*1e-3), '%.3f'%(s[11]*1e-3), '%.1f'%(e[11]*1e-3)],
                ['6V   ','%.2f'%(T[10]-273.15), '%.1f'%(p[10]*1e-3), '%.2f'%(x[10]), '%.1f'%(h[10]*1e-3), '%.3f'%(s[10]*1e-3), '%.1f'%(e[10]*1e-3)],
                ['6IV  ','%.2f'%(T[9]-273.15), '%.1f'%(p[9]*1e-3), '%.2f'%(x[9]), '%.1f'%(h[9]*1e-3), '%.3f'%(s[9]*1e-3), '%.1f'%(e[9]*1e-3)],
                ['6III ','%.2f'%(T[8]-273.15), '%.1f'%(p[8]*1e-3), '%.2f'%(x[8]), '%.1f'%(h[8]*1e-3), '%.3f'%(s[8]*1e-3), '%.1f'%(e[8]*1e-3)],
                ['6II  ','%.2f'%(T[7]-273.15), '%.1f'%(p[7]*1e-3), '%.2f'%(x[7]), '%.1f'%(h[7]*1e-3), '%.3f'%(s[7]*1e-3), '%.1f'%(e[7]*1e-3)],
                ['6I   ','%.2f'%(T[6]-273.15), '%.1f'%(p[6]*1e-3), '%.2f'%(x[6]), '%.1f'%(h[6]*1e-3), '%.3f'%(s[6]*1e-3), '%.1f'%(e[6]*1e-3)],
                ['6    ','%.2f'%(T[5]-273.15), '%.1f'%(p[5]*1e-3), '%.2f'%(x[5]), '%.1f'%(h[5]*1e-3), '%.3f'%(s[5]*1e-3), '%.1f'%(e[5]*1e-3)],
                ['7    ','%.2f'%(T[14]-273.15), '%.1f'%(p[14]*1e-3), '%.2f'%(x[14]), '%.1f'%(h[14]*1e-3), '%.3f'%(s[14]*1e-3), '%.1f'%(e[14]*1e-3)],
                ['7VIII','%.2f'%(T[22]-273.15), '%.1f'%(p[22]*1e-3), '%.2f'%(x[22]), '%.1f'%(h[22]*1e-3), '%.3f'%(s[22]*1e-3), '%.1f'%(e[22]*1e-3)],
                ['7VII ','%.2f'%(T[21]-273.15), '%.1f'%(p[21]*1e-3), '%.2f'%(x[21]), '%.1f'%(h[21]*1e-3), '%.3f'%(s[21]*1e-3), '%.1f'%(e[21]*1e-3)],
                ['7VI  ','%.2f'%(T[20]-273.15), '%.1f'%(p[20]*1e-3), '%.2f'%(x[20]), '%.1f'%(h[20]*1e-3), '%.3f'%(s[20]*1e-3), '%.1f'%(e[20]*1e-3)],
                ['7V   ','%.2f'%(T[19]-273.15), '%.1f'%(p[19]*1e-3), '%.2f'%(x[19]), '%.1f'%(h[19]*1e-3), '%.3f'%(s[19]*1e-3), '%.1f'%(e[19]*1e-3)],
                ['7III ','%.2f'%(T[17]-273.15), '%.1f'%(p[17]*1e-3), '%.2f'%(x[17]), '%.1f'%(h[17]*1e-3), '%.3f'%(s[17]*1e-3), '%.1f'%(e[17]*1e-3)],
                ['7II  ','%.2f'%(T[16]-273.15), '%.1f'%(p[16]*1e-3), '%.2f'%(x[16]), '%.1f'%(h[16]*1e-3), '%.3f'%(s[16]*1e-3), '%.1f'%(e[16]*1e-3)],
                ['7I   ','%.2f'%(T[15]-273.15), '%.1f'%(p[15]*1e-3), '%.2f'%(x[15]), '%.1f'%(h[15]*1e-3), '%.3f'%(s[15]*1e-3), '%.1f'%(e[15]*1e-3)],
                ['9VIII','%.2f'%(T[32]-273.15), '%.1f'%(p[32]*1e-3), '%.2f'%(x[32]), '%.1f'%(h[32]*1e-3), '%.3f'%(s[32]*1e-3), '%.1f'%(e[32]*1e-3)],
                ['9VII ','%.2f'%(T[31]-273.15), '%.1f'%(p[31]*1e-3), '%.2f'%(x[31]), '%.1f'%(h[31]*1e-3), '%.3f'%(s[31]*1e-3), '%.1f'%(e[31]*1e-3)],
                ['9VI  ','%.2f'%(T[30]-273.15), '%.1f'%(p[30]*1e-3), '%.2f'%(x[30]), '%.1f'%(h[30]*1e-3), '%.3f'%(s[30]*1e-3), '%.1f'%(e[30]*1e-3)],
                ['9V   ','%.2f'%(T[29]-273.15), '%.1f'%(p[29]*1e-3), '%.2f'%(x[29]), '%.1f'%(h[29]*1e-3), '%.3f'%(s[29]*1e-3), '%.1f'%(e[29]*1e-3)],
                ['9IV  ','%.2f'%(T[28]-273.15), '%.1f'%(p[28]*1e-3), '%.2f'%(x[28]), '%.1f'%(h[28]*1e-3), '%.3f'%(s[28]*1e-3), '%.1f'%(e[28]*1e-3)],
                ['7IV  ','%.2f'%(T[18]-273.15), '%.1f'%(p[18]*1e-3), '%.2f'%(x[18]), '%.1f'%(h[18]*1e-3), '%.3f'%(s[18]*1e-3), '%.1f'%(e[18]*1e-3)],
                ['9III ','%.2f'%(T[27]-273.15), '%.1f'%(p[27]*1e-3), '%.2f'%(x[27]), '%.1f'%(h[27]*1e-3), '%.3f'%(s[27]*1e-3), '%.1f'%(e[27]*1e-3)],
                ['9II  ','%.2f'%(T[26]-273.15), '%.1f'%(p[26]*1e-3), '%.2f'%(x[26]), '%.1f'%(h[26]*1e-3), '%.3f'%(s[26]*1e-3), '%.1f'%(e[26]*1e-3)],
                ['9I   ','%.2f'%(T[25]-273.15), '%.1f'%(p[25]*1e-3), '%.2f'%(x[25]), '%.1f'%(h[25]*1e-3), '%.3f'%(s[25]*1e-3), '%.1f'%(e[25]*1e-3)],
                ['8    ','%.2f'%(T[23]-273.15), '%.1f'%(p[23]*1e-3), '%.2f'%(x[23]), '%.1f'%(h[23]*1e-3), '%.3f'%(s[23]*1e-3), '%.1f'%(e[23]*1e-3)],
                ['9    ','%.2f'%(T[24]-273.15), '%.1f'%(p[24]*1e-3), '%.2f'%(x[24]), '%.1f'%(h[24]*1e-3), '%.3f'%(s[24]*1e-3), '%.1f'%(e[24]*1e-3)]
              ]
    # print(tabulate(data, headers=header, tablefmt="pretty"))    #tablefmt="latex_booktabs" for latex export
    with open('states.txt', 'w') as file:
        file.write(tabulate(data, headers=header, tablefmt="pretty"))
