import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
plt.rcParams.update({'font.size': 16})

# Capillary pressure parameters
Spc_irr = 0
Spc_max = 1
pc0 = 1e3
alpha = 0.5

Swc = 1e-3
Sgr = 1e-3
S = np.linspace(Swc,1-Sgr,100)

pc = pc0 * np.power(S, -alpha)
dpcdS = -alpha*pc0*np.power(S,-alpha-1)/(Spc_max-Spc_irr)

fig, ax = plt.subplots(1,2, figsize=(10,5))

ax[0].plot(S,pc,c='k',lw=1.5)
ax[0].set_xlabel(r'$S_{w,pc}$ [ - ]')
ax[0].set_ylabel(r'$p_c\left(S_{w,pc}\right)$ [Pa]')
ax[0].grid(True)
ax[1].plot(S,dpcdS,c='k',lw=1.5)
ax[1].set_xlabel(r'$S_{w,pc}$ [ - ]')
ax[1].set_ylabel(r'$\partial p_c/\partial S_{w,pc}$ [Pa]')
# ax[1].set_yscale('log')
ax[1].grid(True)

plt.tight_layout()
plt.savefig('plot_pcFuncs.pdf', dpi=300)
plt.show()