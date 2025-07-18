import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import numpy as np

def pc_model(Sw_pc, pc0, alpha):
    return pc0 * np.power(Sw_pc, -alpha)

plt.style.use('science')
plt.rcParams.update({'font.size': 18})

dt = pd.read_csv('input_pantay.csv')

# fig, ax = plt.subplots(1,2,figsize=(10,5))

# ax[0].scatter(dt['sw'],dt['krw'],c='b',label=r'$k_{rw}$')
# ax[0].scatter(dt['sw'],dt['kro'],c='k',label=r'$k_{ro}$')
# ax[0].set_xlabel(r'$S_w$ [-]')
# ax[0].set_ylabel(r'$k_{ri}$ [-]')
# ax[0].grid(True)
# ax[0].legend()

# ax[1].scatter(dt['sw'],dt['pc'],c='r')
# ax[1].set_xlabel(r'$S_w$ [-]')
# ax[1].set_ylabel(r'$p_c$ [Pa]')
# ax[1].grid(True)

# plt.tight_layout()
# plt.show()

Spc_irr = 0
Spc_max = 1
Sw_pc_exp = (dt['sw'] - Spc_irr) / (Spc_max - Spc_irr)

popt, _ = curve_fit(
    pc_model,
    Sw_pc_exp,
    dt['pc'],
    p0=[1e3, 0.5],            # chute inicial para nw
    bounds=([0, 0.1], [1e7, 10]) 
)
pc0_fit, alpha_fit = popt[0], popt[1]

Sw_fit = np.linspace(np.min(dt['sw']), np.max(dt['sw']), 100)
# Sw_fit = np.linspace(1e-3, 0.999, 100)
Sw_pc_fit = (Sw_fit - Spc_irr) / (Spc_max - Spc_irr)
pc_fit = pc_model(Sw_pc_fit, pc0_fit, alpha_fit)


# Plot para pc
plt.figure(figsize=(5,5))
plt.scatter(Sw_pc_exp, dt['pc'],c='k', label='Panday et al. (1994)')
plt.plot(Sw_pc_fit, pc_fit, label=fr'Fitting $p_c$ ($p_{{c,0}}=${{{pc0_fit:.2f}}}, $\alpha=${{{alpha_fit:.2f}}})', color='red')
plt.xlabel('$S_{w,pc}$')
plt.ylabel('$p_c$')
# plt.yscale('log')
plt.legend(fontsize=12)
plt.grid(True)

plt.tight_layout()
plt.show()









