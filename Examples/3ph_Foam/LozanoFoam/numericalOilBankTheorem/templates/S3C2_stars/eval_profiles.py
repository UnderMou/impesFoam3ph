import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from aux_funcs import *

def S_norm(S,Sr,Swc,Sgr,Sor):
    return np.divide(S-Sr,1-Swc-Sor-Sgr)

import scienceplots
plt.style.use('science')
plt.rcParams.update({'font.size': 22})

residuals = {
    'Swc' : 0.1,
    'Sgr' : 0,
    'Sor' : 0.1
}
Swc = residuals['Swc']
Sgr = residuals['Sgr']
Sor = residuals['Sor']


ti = 0
tf = 60000 
write_interval = 50
t = np.linspace(ti,tf,int(tf/write_interval + 1))
t_find = 53500
idt = np.argmin(np.abs(t-t_find))

# Reading from OpenFOAM
Sas, Sbs, Scs = get_Saturations('.', '.', idt, residuals)

# Reading from FOSSIL
path = 'data_profile_FOSSIL.csv'
data_FOSSIL = pd.read_csv(path)

xD = np.linspace(0,1,len(Sas[0]))
# L = x[-1]
# xD = x/L

# x = data_FOSSIL['Points:0'].to_numpy()
# L = x[-1]
# xD_Fo = x/L

Sg_OF = Sas[0]
So_OF = Scs[0]
# Sg_Fo = data_FOSSIL['Gas_Saturation'].to_numpy()
# So_Fo = data_FOSSIL['Oil_Saturation'].to_numpy()

# Sg_OF = S_norm(Sg_OF,Sgr,Swc,Sgr,Sor)
# So_OF = S_norm(So_OF,Sor,Swc,Sgr,Sor)
# Sg_Fo = S_norm(Sg_Fo,Sgr,Swc,Sgr,Sor)
# So_Fo = S_norm(So_Fo,Sor,Swc,Sgr,Sor)

# Reading Analytical solution from Lyu
data_Lyu_Sg = pd.read_csv('Sg_analytical_LyuC3.2.csv')
data_Lyu_So = pd.read_csv('So_analytical_LyuC3.2.csv')

# Plotting the results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5))

# ax1.plot(xD_Fo, Sg_Fo, linestyle='--', color='b', label='FOSSIL')
ax1.plot(xD, Sg_OF, lw=2,linestyle='--', color='k', label='ImpesFOAM')
ax1.plot(data_Lyu_Sg.iloc[:,0], data_Lyu_Sg.iloc[:,1], color='g', label='Analytical - Tang et al., 2022')
# ax1.set_xlabel(r'$x_D$')
ax1.set_ylabel(r'$S_g$')
# ax1.set_title(r'$S_g$ profile at $T_D=1.0$')
ax1.grid()
ax1.legend(fontsize=18)
# ax1.set_ylim([-0.05,0.85])

# ax2.plot(xD_Fo, So_Fo, linestyle='--', color='b', label='FOSSIL')
ax2.plot(xD, So_OF, lw=2,linestyle='--', color='k', label='ImpesFOAM')
ax2.plot(data_Lyu_So.iloc[:,0], data_Lyu_So.iloc[:,1], color='g')
# ax2.set_xlabel(r'$x_D$')
ax2.set_ylabel(r'$S_o$')
# ax2.set_title(r'$S_o$ profile at $T_D=1.0$')
ax2.grid()
# ax2.set_ylim([-0.05,0.85])

plt.tight_layout()
plt.savefig('LyuC3.2.pdf', dpi=300)
plt.show()
