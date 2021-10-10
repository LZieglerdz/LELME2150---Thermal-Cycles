
#
#===IMPORT PACKAGES============================================================
#

import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from GT_group_42 import gas_turbine, get_stoech, get_lambda, get_lambda2
from ST_simple_group_42 import steam_turbine


#
#===BRAYTON CYCLE==============================================================
#
def brayton():
    P_e = 230e6                     # [W]
    p_1, T_1 = 1e+5, 273.15+15      # [Pa], [K]
    T_3 = 1273.15                   # [K]
    r, eta_pi_c = 18, .90           # [-], [-]
    T_3, k_cc = 273.15+1400, .95    # [K], [-]
    eta_pi_t = .90                  # [-]
    k_mec = .015                    # [-]
    options = p_1, T_1, r, eta_pi_c, T_3, k_cc, eta_pi_t, k_mec
    display = False

    ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG = gas_turbine(P_e,options,display)
    p,T,h,s,e = DAT

    #
    # W = (h[2]-h[3])*eta_mec_t - (h[1]-h[0])/eta_mec_c
    # Q = h[2]-h[1]
    # eta_cyclen = W/Q
    # eta_mec = 1
    # eta_gen = 1
    # eta_en = eta_cyclen*eta_mec*eta_gen
    header_sol  = ['States', 'p\nkPa', 't\n°C', 'h\nkJ/kg', 's\nkJ/kgK', 'e\nkJ/kg']
    data_sol    = [ ['1',100 ,15  ,15.1  ,0.054 ,0    ],
                    ['2',1800,429 ,443.4 ,0.142,403.2 ],
                    ['3',1710,1400,1681.8,1.265,1317.9],
                    ['4',100 ,654 ,731.5 ,1.348,343.9 ] ]
    print(tabulate(data_sol, headers=header_sol, tablefmt="pretty"))    #tablefmt="latex_booktabs" for latex export

    header  = ['States', 'p\nkPa', 't\n°C', 'h\nkJ/kg', 's\nkJ/kgK', 'e\nkJ/kg']
    data    = [ ['1','%.0f'%(p[0]/1000), '%.0f'%(T[0]-273.15), '%.1f'%(h[0]/1000), '%.3f'%(s[0]/1000), '%.0f'%(e[0]/1000)],
                ['2','%.0f'%(p[1]/1000), '%.0f'%(T[1]-273.15), '%.1f'%(h[1]/1000), '%.3f'%(s[1]/1000), '%.0f'%(e[1]/1000)],
                ['3','%.0f'%(p[2]/1000), '%.0f'%(T[2]-273.15), '%.1f'%(h[2]/1000), '%.3f'%(s[2]/1000), '%.0f'%(e[2]/1000)],
                ['4','%.0f'%(p[3]/1000), '%.0f'%(T[3]-273.15), '%.1f'%(h[3]/1000), '%.3f'%(s[3]/1000), '%.0f'%(e[3]/1000)] ]
    print(tabulate(data, headers=header, tablefmt="pretty"))    #tablefmt="latex_booktabs" for latex export

    # print('     | p\t| T \t\t| s-s1\t\t| h-h1 \t\t\t| W : %.2f [kj/kg]' %(W/1000))
    # print('     | [kPa]\t| [°C] \t\t| [kJ/kg.K]\t| [kJ/kg] \t\t| Q : %.2f [kj/kg]' %(Q/1000))
    # print(72*'-'+'|')
    # print('  1  | %.2f\t| %.1f \t\t| %.2f \t\t| %.2f \t\t\t| eta_en : %.3f' %(p[0]/1000, T[0]-273.15, s[0]/1000, h[0]/1000,eta_en))
    # print('  2  | %.2f\t| %.1f \t| %.2f \t\t| %.2f \t\t| eta_pi : %.2f' %(p[1]/1000, T[1]-273.15, s[1]/1000, h[1]/1000,eta_pi))
    # print('  3  | %.2f\t| %.1f \t| %.2f \t\t| %.2f \t\t| eta_mec_c : %.2f' %(p[2]/1000, T[2]-273.15, s[2]/1000, h[2]/1000, eta_mec_c))
    # print('  4  | %.2f\t| %.1f \t| %.2f \t\t| %.2f \t\t| ' %(p[3]/1000, T[3]-273.15, s[3]/1000, h[3]/1000))
    # print(72*'-'+'|')

if __name__ == "__main__":
    brayton()
    print(get_lambda('CH4', 1, 4, 0, 273.15+429, 273.15+1400, 18e5, 17.1e5, 431.5e3,  289.0399801464282, 200, 1))
    print(get_lambda2('CH4', 1, 4, 0, 273.15+429, 273.15+1400, 18e5, 17.1e5, 431.5e3,  289.0399801464282, 200, 1))
    # w = get_steoch('CH4', 1, 4, 0, 273.15+429, 273.15+1400, 18e5, 17.1e5, 289.0399801464282, 200, 1)
    print(17.12043124749108, 2.5170126887435957)
