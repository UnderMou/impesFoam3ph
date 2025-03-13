import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def S_norm(S,Sr,Swc,Sgr,Sor):
    return np.divide(S-Sr,1-Swc-Sor-Sgr)

Swc = 0.1
Sgr = 0
Sor = 0.1

# Reading from OpenFOAM
path = 'data_profile.csv'
data_OpenFOAM = pd.read_csv(path)

# Reading from FOSSIL
path = 'data_profile_FOSSIL.csv'
data_FOSSIL = pd.read_csv(path)

x = data_OpenFOAM['Points:0'].to_numpy()
L = x[-1]
xD = x/L

Sg_OF = data_OpenFOAM['Sa'].to_numpy()
So_OF = data_OpenFOAM['Sc'].to_numpy()
Sg_Fo = data_FOSSIL['Gas_Saturation'].to_numpy()
So_Fo = data_FOSSIL['Oil_Saturation'].to_numpy()

Sg_OF = S_norm(Sg_OF,Sgr,Swc,Sgr,Sor)
So_OF = S_norm(So_OF,Sor,Swc,Sgr,Sor)
Sg_Fo = S_norm(Sg_Fo,Sgr,Swc,Sgr,Sor)
So_Fo = S_norm(So_Fo,Sor,Swc,Sgr,Sor)

# Reading Analytical solution from Lyu
data_Lyu_Sg = pd.read_csv('Sg_analytical_LyuC2.csv')
data_Lyu_So = pd.read_csv('So_analytical_LyuC2.csv')

# Plotting the results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5))

ax1.plot(xD, Sg_Fo, linestyle='--', color='b', label='FOSSIL')
ax1.plot(xD, Sg_OF, linestyle='--', color='k', label='IMPES - OpenFOAM')
ax1.plot(data_Lyu_Sg.iloc[:,0], data_Lyu_Sg.iloc[:,1], color='g', label='Analytical - Lyu, X. et.al., SPE Journal, 2021')
ax1.set_xlabel(r'$x_D$')
ax1.set_ylabel(r'$S_g$')
ax1.set_title(r'$S_g$ profile at $T_D=0.2$')
ax1.grid()
ax1.legend()
# ax1.set_ylim([-0.05,0.85])

ax2.plot(xD, So_Fo, linestyle='--', color='b', label='FOSSIL')
ax2.plot(xD, So_OF, linestyle='--', color='k', label='IMPES - OpenFOAM')
ax2.plot(data_Lyu_So.iloc[:,0], data_Lyu_So.iloc[:,1], color='g')
ax2.set_xlabel(r'$x_D$')
ax2.set_ylabel(r'$S_o$')
ax2.set_title(r'$S_o$ profile at $T_D=0.2$')
ax2.grid()
# ax2.set_ylim([-0.05,0.85])

plt.tight_layout()
plt.savefig('LyuC2.png', dpi=300)
plt.show()
