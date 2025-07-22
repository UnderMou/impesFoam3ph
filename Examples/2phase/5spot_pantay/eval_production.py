import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import numpy as np

plt.style.use('science')

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 18

data = pd.read_csv("production.txt", delimiter=',')
print(data.head())
print(data.tail())

# Production Pantay data
pantay_experiment = pd.read_csv('production_pantay_experimental.csv', header=None)
pantay_numericalFe = pd.read_csv('production_pantay_numericalFE.csv', header=None)

# Domain info
L = 142.245 # [m]
h = 10 # [m]
V = L**2 * h # [m3]
phi = 0.2 # [-]
PV = V*phi

t = data['TimeStep'].to_numpy()
dt = t[1:] - t[:-1]
qo_prod = data['qo_prod'].to_numpy() # [m3/s]  
qw_inj = data['qw_inj'].to_numpy() # [m3/s]

Qw_inj = np.zeros_like(qw_inj)
Qo_prod = np.zeros_like(qo_prod)

for i in range(len(Qw_inj)):
    Qw_inj[i] = np.sum(np.multiply(dt[:i],qw_inj[:i]))  # [m3] 
    Qo_prod[i] = np.sum(np.multiply(dt[:i],qo_prod[:i])) # [m3]

# plt.figure(figsize=(7,7))
# plt.plot(t,qw,c='b',label=r"$Q_w$")
# plt.plot(t,qo,c='k',label=r"$Q_o$")
# plt.grid()
# plt.legend()
# # plt.xlim([0,t[-1]+3000])
# # plt.xlim([0,300000])
# # plt.ylim([0,4])
# plt.xlabel('time [s]')
# plt.ylabel('production []')
# plt.show()

# INJECTED WATER

plt.figure(figsize=(10,5))
plt.plot(Qw_inj/PV,Qo_prod/PV,c='k',label=r"Numerical - ImpesFOAM")
# plt.plot(t,Qo_prod/PV,c='k',label=r"Numerical - ImpesFOAM")

plt.plot(pantay_experiment.iloc[:,0], pantay_experiment.iloc[:,1], c='b', label='Experiment - Gaucher and Lindley (1960)')
plt.plot(pantay_numericalFe.iloc[:,0], pantay_numericalFe.iloc[:,1], c='r', label='Numerical - FE - Panday et al. (1994)')
plt.grid()
plt.legend(loc='lower right')
plt.xlabel(r'Water injected [PV]')
plt.ylabel(r'Cumulative oil recovery [PV]')
plt.xlim([-0.1,2.0])

plt.tight_layout()
plt.savefig('production_5spots.pdf', dpi=300)
plt.show()