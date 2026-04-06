import uqtopus as uqt
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import ternary
import matplotlib 
import pandas as pd

exp_dir = '.'

t = np.arange(0.1, 20, 0.1)
nx = 1000


data_single = uqt.parse_openfoam_case(
    case_dir=exp_dir,
    variables=["Sa", "Sb", "Sc", "Cs", "MRF"],
    time_dirs=[f"{i:g}" for i in t]
)

L = 30.48
x = np.linspace(0, L, nx)
xD = x/L
u = 1.0 # m/s
tD = u*t/L
tD_idx = np.argmin(np.abs(tD - 0.1))

Sc = data_single.Sc[tD_idx,:]
Sb = data_single.Sb[tD_idx,:]
Sa = data_single.Sa[tD_idx,:]
Cs = data_single.Cs[tD_idx,:]
MRF = data_single.MRF[tD_idx,:]

vD = xD/tD[tD_idx]
    
plt.plot(xD, Sa, label="Sa", color="green")
plt.ylim([0,1])
plt.show()

datafr = {'So': Sc,
          'Sw': Sb,
          'Sg': Sa,
          'Cs': Cs,
          'MRF': MRF,
          'xD': xD}
datafr = pd.DataFrame(datafr)

filename = f'profiles_secondaryBaseCase_tD_{tD[tD_idx]:.2f}.csv'
datafr.to_csv(os.path.join(exp_dir, filename), index=False)

X = pd.read_csv(os.path.join(exp_dir, filename))
print(X.head())

plt.plot(X.xD, X.Cs, label="Cs", color="red")
plt.plot(X.xD, X.So, label="So", color="black")
plt.plot(X.xD, X.Sw, label="Sw", color="blue")
plt.plot(X.xD, X.Sg, label="Sg", color="green")
plt.show()
plt.close()
plt.plot(X.xD, X.MRF, label="MRF", color="red")
plt.show()