import CoolProp
import numpy as np


from CoolProp.Plots.SimpleCycles import StateContainer

T0 = 300.000; p0 = 200000.000; h0 = 112745.749; s0 = 393.035

cycle_states = StateContainer()
#
cycle_states[0,'H'] = h0

cycle_states[0]['S'] = s0

cycle_states[0][CoolProp.iP] = p0

cycle_states[0,CoolProp.iT] = T0
#
# cycle_states[1,"T"] = 300.064
# cycle_states[2,"T"] = 300.064

s_points = np.array( ( 0,1,2,3,4 ))#*1e-3
h_points = np.array( ( 0,1,2,3,4 ))#*1e-3
T_points = np.array( ( 0,1,2,3,4 ))
for i in range(len(s_points)):
    cycle_states[i,'S'] = s_points[i]
    cycle_states[i,'H'] = h_points[i]
    cycle_states[i,'T'] = T_points[i]
print(cycle_states)
