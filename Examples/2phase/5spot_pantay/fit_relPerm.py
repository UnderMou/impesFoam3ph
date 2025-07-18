import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import numpy as np

def krw_model(Seff, nw):
    return 0.562 * Seff**nw

def kro_model(Seff, no):
    return 0.813 * (1-Seff)**no

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

Swc = 0.3
Sor = 0.067
Seff_exp = (dt['sw'] - Swc) / (1 - Swc - Sor)

popt, _ = curve_fit(
    krw_model,
    Seff_exp,
    dt['krw'],
    p0=[.5],            # chute inicial para nw
    bounds=([.1], [10.0])  # limites para nw
)
nw_fit = popt[0]

popt, _ = curve_fit(
    kro_model,
    Seff_exp,
    dt['kro'],
    p0=[.5],            # chute inicial para nw
    bounds=([.1], [10.0])  # limites para nw
)
no_fit = popt[0]



Sw_fit = np.linspace(np.min(dt['sw']), np.max(dt['sw']), 100)
Seff_fit = (Sw_fit - Swc) / (1 - Swc - Sor)
krw_fit = krw_model(Seff_fit, nw_fit)
kro_fit = kro_model(Seff_fit, no_fit)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))

# Plot para krw
ax[0].scatter(Seff_exp, dt['krw'],c='k', label='Panday et al. (1994)')
ax[0].plot(Seff_fit, krw_fit, label=fr'Fitting $k_{{rw}}$ ($n_w=${{ {nw_fit:.2f} }})', color='red')
# ax[0].set_xscale('log')
# ax[0].set_yscale('log')
ax[0].set_xlabel('$S_{eff}$')
ax[0].set_ylabel('$k_{rw}$')
ax[0].legend(fontsize=12)
ax[0].grid(True)
# ax[0].set_title('Ajuste de $k_{rw}$')

# Plot para kro
ax[1].scatter(Seff_exp, dt['kro'],c='k', label='Panday et al. (1994)')
# Se você tiver kro_fit, plote abaixo. Caso contrário, remova essa linha.
ax[1].plot(Seff_fit, kro_fit, label=fr'Fitting $k_{{ro}}$ ($n_o$={{ {no_fit:.2f} }})', color='red')
# ax[1].set_xscale('log')
# ax[1].set_yscale('log')
ax[1].set_xlabel('$S_{eff}$')
ax[1].set_ylabel('$k_{ro}$')
ax[1].legend(fontsize=12)
ax[1].grid(True)
# ax[1].set_title('Ajuste de $k_{ro}$')

plt.tight_layout()
plt.show()









