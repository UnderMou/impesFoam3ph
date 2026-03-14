import uqtopus as uqt
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import ternary
import matplotlib 

exp_dir = '.'

t = np.arange(0.1, 20, 0.1)


data_single = uqt.parse_openfoam_case(
    case_dir=exp_dir,
    variables=["Sa", "Sb", "Sc"],
    time_dirs=[f"{i:g}" for i in t]
)

t_des = 2.5
t_idx = np.argmin(np.abs(np.asarray(data_single.time) - t_des))

Sc = data_single.Sc[t_idx,:]
Sb = data_single.Sb[t_idx,:]
Sa = data_single.Sa[t_idx,:]

L = 30.48
x = np.linspace(0, L, len(Sa))
xD = x/L

u = 1.0 # m/s
tD = u*t/L

tD_idx = np.argmin(np.abs(tD - 0.3))

vD = xD/tD[tD_idx]
    
plt.plot(vD, Sa, label="Sa", color="green")
plt.ylim([0,1])
plt.show()
