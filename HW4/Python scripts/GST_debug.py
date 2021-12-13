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

path = 'GST_group_' + '42';
print('\n '+path)
GST = importlib.import_module(path).GST
P_eg = 225e+6
P_es = 140e+6
comb = 0, 4
options = comb, 288.15, 1e+5, 15+273.15, 100e3, 1250+273.15, 15, 0.9, 0.9, 1-0.05, 0.015, 50+273.15, 11e6, 2.8e6, 565+273.15, 0.4e6, 318+273.15,
6e3, 0.95
display = False
(ETA,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG) = GST(P_eg,P_es,options,display)
# (ETA,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION) = GST(P_eg,P_es,options,display)
