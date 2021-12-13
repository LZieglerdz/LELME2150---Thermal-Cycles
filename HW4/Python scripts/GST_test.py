"""
LELME2150 - Thermal cycles
Homework 4 - GST test code

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
for gr in GROUPS:
    try:
        path = 'GST_group_' + gr;
        print('\n '+path)
        GST = importlib.import_module(path).GST
        P_eg = 225e+6
        P_es = 140e+6
        comb = 0, 4
        options = comb, 288.15, 1e+5, 15+273.15, 100e3, 1250+273.15, 15, 0.9, 0.9, 1-0.05, 0.015, 50+273.15, 11e6, 2.8e6, 565+273.15, 0.4e6, 318+273.15, 6e3, 0.95
        display = False
        try:
            (ETA,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG) = GST(P_eg,P_es,options,display)
            # (ETA,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION) = GST(P_eg,P_es,options,display)
            (LHV,e_c,excess_air,cp_gas,gas_prop) = COMBUSTION
            (dotm_a,dotm_f,dotm_g,dotm_v, dotm_vLP, dotm_vIP, dotm_vHP) = MASSFLOW
            # (eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex) = ETA
            (loss_mec,loss_cond,loss_chimney) = DATEN
            (loss_mec,loss_rotex,loss_combex,loss_chemex,loss_transex,loss_totex,loss_condex) = DATEX
            (fig_pie_en,fig_pie_ex,fig_Ts_diagram,fig_hs_diagram,fig_heat_exchange) = FIG
            (p,T,h,s,e,x) = DAT
            (p_1g,p_2g,p_3g,p_4g,p_5g,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_8p,p_8pp,p_9,p_9p,p_9pp,p_10p,p_10pp) = p
            (T_1g,T_2g,T_3g,T_4g,T_5g,T_1,T_2,T_3,T_4,T_5,T_6,T_7,T_8,T_8p,T_8pp,T_9,T_9p,T_9pp,T_10p,T_10pp) = T
            (s_1g,s_2g,s_3g,s_4g,s_5g,s_1,s_2,s_3,s_4,s_5,s_6,s_7,s_8,s_8p,s_8pp,s_9,s_9p,s_9pp,s_10p,s_10pp) = s
            (h_1g,h_2g,h_3g,h_4g,h_5g,h_1,h_2,h_3,h_4,h_5,h_6,h_7,h_8,h_8p,h_8pp,h_9,h_9p,h_9pp,h_10p,h_10pp) = h
            (e_1g,e_2g,e_3g,e_4g,e_5g,e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_8p,e_8pp,e_9,e_9p,e_9pp,e_10p,e_10pp) = e
            (x_1g,x_2g,x_3g,x_4g,x_5g,x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_8p,x_8pp,x_9,x_9p,x_9pp,x_10p,x_10pp) = x

            print(50*'-')
            print('COMBUSTION: ', COMBUSTION)
            print(50*'-')
            print('MASSFLOW: ', MASSFLOW)
            print(50*'-')
            print('ETA: ', ETA)
            print(50*'-')
            print('DATEN: ',DATEN)
            print(50*'-')
            print('DATEX: ', DATEX)
            print(50*'-')

            header  = ['States', 't\n°C', 'p\nkPa', 'x\n-', 'h\nkJ/kg', 's\nkJ/kgK', 'e\nkJ/kg']
            data    = [ ['1g    ','%.2f'%(T[0]-273.15),  '%.1f'%(p[0]*1e-3),  '%.2f'%(x[0]),  '%.1f'%(h[0]*1e-3),  '%.3f'%(s[0]*1e-3),  '%.1f'%(e[0]*1e-3)],
            ['2g    ','%.2f'%(T[1]-273.15),  '%.1f'%(p[1]*1e-3),  '%.2f'%(x[1]),  '%.1f'%(h[1]*1e-3),  '%.3f'%(s[1]*1e-3),  '%.1f'%(e[1]*1e-3)],
            ['3g    ','%.2f'%(T[2]-273.15),  '%.1f'%(p[2]*1e-3),  '%.2f'%(x[2]),  '%.1f'%(h[2]*1e-3),  '%.3f'%(s[2]*1e-3),  '%.1f'%(e[2]*1e-3)],
            ['4g    ','%.2f'%(T[3]-273.15),  '%.1f'%(p[3]*1e-3),  '%.2f'%(x[3]),  '%.1f'%(h[3]*1e-3),  '%.3f'%(s[3]*1e-3),  '%.1f'%(e[3]*1e-3)],
            ['5g    ','%.2f'%(T[4 ]-273.15), '%.1f'%(p[4 ]*1e-3), '%.2f'%(x[4 ]), '%.1f'%(h[4 ]*1e-3), '%.3f'%(s[4 ]*1e-3), '%.1f'%(e[4 ]*1e-3)],
            ['1     ','%.2f'%(T[5 ]-273.15), '%.1f'%(p[5 ]*1e-3), '%.2f'%(x[5 ]), '%.1f'%(h[5 ]*1e-3), '%.3f'%(s[5 ]*1e-3), '%.1f'%(e[5 ]*1e-3)],
            ['2     ','%.2f'%(T[6 ]-273.15), '%.1f'%(p[6 ]*1e-3), '%.2f'%(x[6 ]), '%.1f'%(h[6 ]*1e-3), '%.3f'%(s[6 ]*1e-3), '%.1f'%(e[6 ]*1e-3)],
            ['3     ','%.2f'%(T[7 ]-273.15), '%.1f'%(p[7 ]*1e-3), '%.2f'%(x[7 ]), '%.1f'%(h[7 ]*1e-3), '%.3f'%(s[7 ]*1e-3), '%.1f'%(e[7 ]*1e-3)],
            ['4     ','%.2f'%(T[8 ]-273.15), '%.1f'%(p[8 ]*1e-3), '%.2f'%(x[8 ]), '%.1f'%(h[8 ]*1e-3), '%.3f'%(s[8 ]*1e-3), '%.1f'%(e[8 ]*1e-3)],
            ['5     ','%.2f'%(T[9 ]-273.15), '%.1f'%(p[9 ]*1e-3), '%.2f'%(x[9 ]), '%.1f'%(h[9 ]*1e-3), '%.3f'%(s[9 ]*1e-3), '%.1f'%(e[9 ]*1e-3)],
            ['6     ','%.2f'%(T[10]-273.15), '%.1f'%(p[10]*1e-3), '%.2f'%(x[10]), '%.1f'%(h[10]*1e-3), '%.3f'%(s[10]*1e-3), '%.1f'%(e[10]*1e-3)],
            ['7     ','%.2f'%(T[11]-273.15), '%.1f'%(p[11]*1e-3), '%.2f'%(x[11]), '%.1f'%(h[11]*1e-3), '%.3f'%(s[11]*1e-3), '%.1f'%(e[11]*1e-3)],
            ['8     ','%.2f'%(T[12]-273.15), '%.1f'%(p[12]*1e-3), '%.2f'%(x[12]), '%.1f'%(h[12]*1e-3), '%.3f'%(s[12]*1e-3), '%.1f'%(e[12]*1e-3)],
            ['8\'   ','%.2f'%(T[13]-273.15), '%.1f'%(p[13]*1e-3), '%.2f'%(x[13]), '%.1f'%(h[13]*1e-3), '%.3f'%(s[13]*1e-3), '%.1f'%(e[13]*1e-3)],
            ['8\"   ','%.2f'%(T[14]-273.15), '%.1f'%(p[14]*1e-3), '%.2f'%(x[14]), '%.1f'%(h[14]*1e-3), '%.3f'%(s[14]*1e-3), '%.1f'%(e[14]*1e-3)],
            ['9     ','%.2f'%(T[15]-273.15), '%.1f'%(p[15]*1e-3), '%.2f'%(x[15]), '%.1f'%(h[15]*1e-3), '%.3f'%(s[15]*1e-3), '%.1f'%(e[15]*1e-3)],
            ['9\'   ','%.2f'%(T[16]-273.15), '%.1f'%(p[16]*1e-3), '%.2f'%(x[16]), '%.1f'%(h[16]*1e-3), '%.3f'%(s[16]*1e-3), '%.1f'%(e[16]*1e-3)],
            ['9\"   ','%.2f'%(T[17]-273.15), '%.1f'%(p[17]*1e-3), '%.2f'%(x[17]), '%.1f'%(h[17]*1e-3), '%.3f'%(s[17]*1e-3), '%.1f'%(e[17]*1e-3)],
            ['10\'  ','%.2f'%(T[18]-273.15), '%.1f'%(p[18]*1e-3), '%.2f'%(x[18]), '%.1f'%(h[18]*1e-3), '%.3f'%(s[18]*1e-3), '%.1f'%(e[18]*1e-3)],
            ['10\"  ','%.2f'%(T[19]-273.15), '%.1f'%(p[19]*1e-3), '%.2f'%(x[19]), '%.1f'%(h[19]*1e-3), '%.3f'%(s[19]*1e-3), '%.1f'%(e[19]*1e-3)]
            ]
            # print(tabulate(data, headers=header, tablefmt="pretty"))    #tablefmt="latex_booktabs" for latex export
            with open('states.txt', 'w') as file:
                file.write(tabulate(data, headers=header, tablefmt="pretty"))

            print('\n Basic ST passed: your code will be graded!')
        except:
            print('\n Basic ST failed: please correct your code!')
    except:
        print('\n Basic ST failed: please correct your code!')
