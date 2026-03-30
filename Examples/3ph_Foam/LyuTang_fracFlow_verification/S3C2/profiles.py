import uqtopus as uqt
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import ternary
import matplotlib 
import pandas as pd
import scienceplots

def S_norm(S,Sr,Swc,Sgr,Sor):
    return np.divide(S-Sr,1-Swc-Sor-Sgr)

lozanoPath = "./"

sg_lozano = pd.read_csv(os.path.join(lozanoPath, "Sg_analytical_LyuC3.2.csv"), header=0)
so_lozano = pd.read_csv(os.path.join(lozanoPath, "So_analytical_LyuC3.2.csv"), header=0)

x_interp = np.linspace(0,1,1000)
sg_lozano_interp = np.interp(x_interp, sg_lozano.iloc[:,0], sg_lozano.iloc[:,1])
so_lozano_interp = np.interp(x_interp, so_lozano.iloc[:,0], so_lozano.iloc[:,1])
sw_lozano_interp = 1 - sg_lozano_interp - so_lozano_interp

exp_dir = '.'

t = np.arange(50, 60000, 50)


data_single = uqt.parse_openfoam_case(
    case_dir=exp_dir,
    variables=["Sa", "Sb", "Sc"],
    time_dirs=[f"{i:g}" for i in t]
)

t_des = 53500
t_idx = np.argmin(np.abs(np.asarray(data_single.time) - t_des))

Sc = data_single.Sc[t_idx,:]
Sb = data_single.Sb[t_idx,:]
Sa = data_single.Sa[t_idx,:]

L = 1
x = np.linspace(0, L, len(Sa))
xD = x/L

u = 1.4111e-5 # m/s
tD = u*t/L

tD_idx = np.argmin(np.abs(tD - 1.0))

vD = xD/tD[tD_idx]
Sc = data_single.Sc[t_idx,:]
Sb = data_single.Sb[t_idx,:]
Sa = data_single.Sa[t_idx,:]


Swc = 0.1
Sgr = 0
Sor = 0.1
Sa = S_norm(Sa,Sgr,Swc,Sgr,Sor)
Sb = S_norm(Sb,Swc,Swc,Sgr,Sor)
Sc = S_norm(Sc,Sor,Swc,Sgr,Sor)
    
# plt.scatter(xD, Sa, label="Sa", color="green", s=10)
# plt.scatter(xD, Sb, label="Sb", color="blue", s=10)
# plt.scatter(xD, Sc, label="Sc", color="black", s=10)

# plt.plot(sg_lozano.iloc[:,0], sg_lozano.iloc[:,1], label="Lozano SG", color="green")
# plt.plot(sw_lozano.iloc[:,0], sw_lozano.iloc[:,1], label="Lozano SW", color="blue")
# plt.plot(so_lozano.iloc[:,0], so_lozano.iloc[:,1], label="Lozano SO", color="black")
# plt.ylim([0,1])
# plt.show()

plt.style.use(['science', 'grid'])   # science style + subtle grid

matplotlib.rcParams.update({'font.size': 22})


fig, ax = plt.subplots(figsize=(8,9))


# Reference (Lozano)
ax.plot(x_interp, sg_lozano_interp,
        label="$S_g$, Lyu et al. 2021", color="green", lw=6, alpha=0.7)

ax.plot(x_interp, sw_lozano_interp,
        label="$S_w$, Lyu et al. 2021", color="blue", lw=6, alpha=0.7)

ax.plot(x_interp, so_lozano_interp,
        label="$S_o$, Lyu et al. 2021", color="black", lw=6, alpha=0.7)

# Simulation data
ax.plot(xD, Sa, label="$S_g$, Numerical", color="lime", lw=4, linestyle='dashed')
ax.plot(xD, Sb, label="$S_w$, Numerical", color="cornflowerblue",  lw=4, linestyle='dashed')
ax.plot(xD, Sc, label="$S_o$, Numerical", color="dimgray", lw=4, linestyle='dashed')

# Formatting
ax.set_ylim(0,1)
ax.set_xlabel(r"$x_D$")
ax.set_ylabel("Saturation")

# Legend outside (optimized layout)
ax.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, 1.35),
    ncol=2,
    frameon=False,
    columnspacing=1.2,
    handlelength=2.5,
    handletextpad=0.5
)

ax.autoscale(tight=True)

plt.tight_layout()
plt.savefig("profilesLyu.png", dpi=300)  # Save the figure with high resolution
plt.show()
