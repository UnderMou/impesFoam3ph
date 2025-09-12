import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import numpy as np
from matplotlib.lines import Line2D

plt.style.use('science')

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 22

# Domain info
Soi = 0.28   # Original oil in place OOIP
Lx = 360     # [m]
Ly = 100     # [m]
h = 30          # [m]
V = Lx*Ly * h   # [m3]
phi = 0.2       # [-]
PV = V*phi      # [m3]

# PLOT
paths = ['noFoam_G', 'withFoam_G']
linestyles = ['--', '-']
labels = ['without foam', 'with foam']
fig, ax = plt.subplots(1,3,figsize=(20,6))

for j in range(len(paths)):
    data = pd.read_csv(paths[j] + "/production.txt", delimiter=',')
    print(data.head())
    print(data.tail())

    t = data['TimeStep'].to_numpy()
    print('t_final: %.2f years' % (t[-1]/(60*60*24*365)))

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

    # CHECK REPORT
    Qt_inj = Qw_inj[-1]+Qg_inj[-1]
    print('Total injected: ', Qt_inj)
    Qt_prod = Qw_prod[-1]+Qg_prod[-1]+Qo_prod[-1]
    print('Total produced: ', Qt_prod)
    print('Cons.mass error:', 100*abs(Qt_inj-Qt_prod)/Qt_inj, '%')

    # ORF
    OOIP = Soi
    ORF = ((Qo_prod/PV) / OOIP) * 100

    # ax[0].plot((Qw_inj+Qg_inj)/PV,Qo_prod/PV,linestyle = linestyles[j], c='k',lw=2,label='oil '+labels[j])
    # ax[0].plot((Qw_inj+Qg_inj)/PV,Qw_prod/PV,linestyle = linestyles[j], c='b',lw=2,label='water '+labels[j])
    # ax[0].plot((Qw_inj+Qg_inj)/PV,Qg_prod/PV,linestyle = linestyles[j], c='g',lw=2,label='gas '+labels[j])
    ax[0].plot((Qw_inj+Qg_inj)/PV, Qo_prod/PV,
               linestyle=linestyles[j], c='k', lw=2)
    ax[0].plot((Qw_inj+Qg_inj)/PV, Qw_prod/PV,
               linestyle=linestyles[j], c='b', lw=2)
    ax[0].plot((Qw_inj+Qg_inj)/PV, Qg_prod/PV,
               linestyle=linestyles[j], c='g', lw=2)
    
    ax[1].plot((Qw_inj+Qg_inj)/PV,Qg_prod/Qo_prod,linestyle = linestyles[j],c='k',lw=2,label=labels[j])
      
    ax[2].plot((Qw_inj+Qg_inj)/PV,ORF,linestyle = linestyles[j],c='k',lw=2,label=labels[j])
    
ax[0].grid()
ax[0].set_xlabel(r'Water+Gas injected [PV]')
ax[0].set_ylabel(r'Cumulative phase recovery [PV]')
# ax[0].set_ylim([-0.02,0.7])
# ax[0].set_xlim([-0.25,1.35])
# ax[0].legend()
# Legend for fluids
fluid_legend = [
    Line2D([0], [0], color='k', lw=2, label='Oil'),
    Line2D([0], [0], color='b', lw=2, label='Water'),
    Line2D([0], [0], color='g', lw=2, label='Gas'),
]
legend1 = ax[0].legend(handles=fluid_legend, loc="upper left")
ax[0].add_artist(legend1)

# Legend for foam
foam_legend = [
    Line2D([0], [0], color='k', lw=2, linestyle='--', label='Without foam'),
    Line2D([0], [0], color='k', lw=2, linestyle='-',  label='With foam'),
    Line2D([0], [0], color='m', lw=2, linestyle='dotted', label='GOR $=1$')
]
ax[1].legend(handles=foam_legend, loc="upper left")

ax[1].set_xlabel(r'Water+Gas injected [PV]')
ax[1].set_ylabel(r'GOR [-]')
ax[1].axhline(1,linestyle='dotted',c='m',lw=2,label='GOR $=1$')
ax[1].grid()
# ax[1].legend()

ax[2].grid()
ax[2].set_xlabel(r'Water+Gas injected [PV]')
ax[2].set_ylabel(r'ORF [\%]')
# ax[2].legend()

plt.tight_layout()
plt.savefig('production_spe.pdf', dpi=300)
plt.show()