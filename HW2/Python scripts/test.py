
#
#===IMPORT PACKAGES============================================================
#

import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from GT_group_42 import gas_turbine, get_lambda
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

    header_sol  = ['States', 'p\nkPa', 't\n°C', 'h\nkJ/kg', 's\nkJ/kgK', 'e\nkJ/kg']
    data_sol    = [ ['1',100 ,15  ,15.1  ,0.054 ,0    ],
                    ['2',1800,429 ,443.4 ,0.142,403.2 ],
                    ['3',1710,1400,1681.8,1.265,1317.9],
                    ['4',100 ,654 ,731.5 ,1.348,343.9 ] ]
    # print(tabulate(data_sol, headers=header_sol, tablefmt="pretty"))    #tablefmt="latex_booktabs" for latex export

    header  = ['States', 'p\nkPa', 't\n°C', 'h\nkJ/kg', 's\nkJ/kgK', 'e\nkJ/kg']
    data    = [ ['1','%.0f'%(p[0]/1000), '%.0f'%(T[0]-273.15), '%.1f'%(h[0]/1000), '%.3f'%(s[0]/1000), '%.0f'%(e[0]/1000)],
                ['2','%.0f'%(p[1]/1000), '%.0f'%(T[1]-273.15), '%.1f'%(h[1]/1000), '%.3f'%(s[1]/1000), '%.0f'%(e[1]/1000)],
                ['3','%.0f'%(p[2]/1000), '%.0f'%(T[2]-273.15), '%.1f'%(h[2]/1000), '%.3f'%(s[2]/1000), '%.0f'%(e[2]/1000)],
                ['4','%.0f'%(p[3]/1000), '%.0f'%(T[3]-273.15), '%.1f'%(h[3]/1000), '%.3f'%(s[3]/1000), '%.0f'%(e[3]/1000)] ]
    print(tabulate(data, headers=header, tablefmt="pretty"))    #tablefmt="latex_booktabs" for latex export

    # print(ETA)
    # print('eta_cyclen: %.2f,\neta_toten %.2f,\neta_cyclex %.2f,\neta_totex %.2f,\neta_rotex %.2f,\neta_combex %.2f' %ETA)

if __name__ == "__main__":
    brayton()
