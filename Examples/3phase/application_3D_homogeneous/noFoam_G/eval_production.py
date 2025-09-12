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

# # Production Pantay data
# pantay_experiment = pd.read_csv('production_pantay_experimental.csv', header=None)
# pantay_numericalFe = pd.read_csv('production_pantay_numericalFE.csv', header=None)
# pantay_numericalFd = pd.read_csv('production_pantay_numericalFD.csv', header=None)

# Domain info
Lx = 360     # [m]
Ly = 100     # [m]
Lz = 85      # [m]
V = Lx*Ly*Lz    # [m3]
phi = 0.2       # [-]
PV = V*phi      # [m3]

t = data['TimeStep'].to_numpy()
dt = t[1:] - t[:-1]
qo_prod = data['qo_prod'].to_numpy() # [m3/s]  
qw_inj = data['qw_inj'].to_numpy()   # [m3/s]
qw_prod = data['qw_prod'].to_numpy() # [m3/s]
qg_inj = data['qg_inj'].to_numpy()   # [m3/s]
qg_prod = data['qg_prod'].to_numpy() # [m3/s]

Qw_inj = np.zeros_like(qw_inj)
Qw_prod = np.zeros_like(qw_inj)
Qo_prod = np.zeros_like(qo_prod)
Qg_inj = np.zeros_like(qg_inj)
Qg_prod = np.zeros_like(qg_inj)

# CUMULATIVE PRODUCTION
for i in range(len(Qw_inj)):
    Qw_inj[i] = np.sum(np.multiply(dt[:i],qw_inj[:i]))      # [m3] 
    Qo_prod[i] = np.sum(np.multiply(dt[:i],qo_prod[:i]))    # [m3]
    Qw_prod[i] = np.sum(np.multiply(dt[:i],qw_prod[:i]))    # [m3] 
    Qg_inj[i] = np.sum(np.multiply(dt[:i],qg_inj[:i]))      # [m3] 
    Qg_prod[i] = np.sum(np.multiply(dt[:i],qg_prod[:i]))

# ORF
Soi = 0.28
OOIP = Soi
ORF = ((Qo_prod/PV) / OOIP) * 100


# PLOT
fig, ax = plt.subplots(1,3,figsize=(10,5))

# plt.plot(pantay_experiment.iloc[:,0], pantay_experiment.iloc[:,1], c='b', label='Experiment - Gaucher and Lindley (1960)',lw=2)
# plt.plot(pantay_numericalFe.iloc[:,0], pantay_numericalFe.iloc[:,1], c='r', label='Numerical - FE - Panday et al. (1994)',lw=2)
# plt.plot(pantay_numericalFd.iloc[:,0], pantay_numericalFd.iloc[:,1], c='g', label='Numerical - FD - Panday et al. (1994)',lw=2)

ax[0].plot((Qw_inj+Qg_inj)/PV,Qo_prod/PV,c='k',lw=2,label='oil')
ax[0].plot((Qw_inj+Qg_inj)/PV,Qw_prod/PV,c='b',lw=2,label='water')
ax[0].plot((Qw_inj+Qg_inj)/PV,Qg_prod/PV,c='g',lw=2,label='gas')
ax[0].grid()
ax[0].set_xlabel(r'Water+Gas injected [PV]')
ax[0].set_ylabel(r'Cumulative phase recovery [PV]')
ax[0].legend()

ax[1].plot((Qw_inj+Qg_inj)/PV,Qg_prod/Qo_prod,c='k',lw=2)
ax[1].grid()
ax[1].set_xlabel(r'Water+Gas injected [PV]')
ax[1].set_ylabel(r'GOR [-]')

ax[2].plot((Qw_inj+Qg_inj)/PV,ORF,c='k',lw=2)
ax[2].grid()
ax[2].set_xlabel(r'Water+Gas injected [PV]')
ax[2].set_ylabel(r'ORF [\%]')

plt.tight_layout()
plt.savefig('production_5spots.pdf', dpi=300)
plt.show()