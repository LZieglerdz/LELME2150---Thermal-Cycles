
#
#===IMPORT PACKAGES============================================================
#

import numpy as np
import matplotlib.pyplot as plt
from basic_cycles_group_42 import gas_turbine
from basic_cycles_group_42 import steam_turbine


#
#===BRAYTON CYCLE==============================================================
#
def brayton():
    # p_1, T_1 = 1e+5, 293.15 # [Pa], [K]
    # p_2 = 1e+7 # [Pa]
    # p_3, T_3 = 1e+7, 1273.15 # [Pa], [K]

    p_1, T_1 = 1e+5, 293.15 # [Pa], [K]
    p_2 = 17.8e+5 # [Pa]
    p_3, T_3 = 17.8e+5, 1273.15 # [Pa], [K]

    eta_pi = 0.90 # [-]
    eta_mec_c, eta_mec_t = 0.98, 0.98 # [-], [-]

    p,T,s,h,eta_en = gas_turbine(p_1,T_1,p_2,p_3,T_3,eta_pi,eta_mec_c,eta_mec_t)
    W = (h[2]-h[3])*eta_mec_t - (h[1]-h[0])/eta_mec_c
    Q = h[2]-h[1]
    eta_cyclen = W/Q
    eta_mec = 1
    eta_gen = 1
    eta_en = eta_cyclen*eta_mec*eta_gen

    print('     | p\t| T \t\t| s-s1\t\t| h-h1 \t\t\t| W : %.2f [kj/kg]' %(W/1000))
    print('     | [kPa]\t| [°C] \t\t| [kJ/kg.K]\t| [kJ/kg] \t\t| Q : %.2f [kj/kg]' %(Q/1000))
    print(72*'-'+'|')
    print('  1  | %.2f\t| %.1f \t\t| %.2f \t\t| %.2f \t\t\t| eta_en : %.3f' %(p[0]/1000, T[0]-273.15, s[0]/1000, h[0]/1000,eta_en))
    print('  2  | %.2f\t| %.1f \t| %.2f \t\t| %.2f \t\t| eta_pi : %.2f' %(p[1]/1000, T[1]-273.15, s[1]/1000, h[1]/1000,eta_pi))
    print('  3  | %.2f\t| %.1f \t| %.2f \t\t| %.2f \t\t| eta_mec_c : %.2f' %(p[2]/1000, T[2]-273.15, s[2]/1000, h[2]/1000, eta_mec_c))
    print('  4  | %.2f\t| %.1f \t| %.2f \t\t| %.2f \t\t| ' %(p[3]/1000, T[3]-273.15, s[3]/1000, h[3]/1000))
    print(72*'-'+'|')

if __name__ == "__main__":
    brayton()
