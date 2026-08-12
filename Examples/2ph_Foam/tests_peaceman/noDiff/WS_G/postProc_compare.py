import uqtopus as uqt
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import ternary
import matplotlib
import pandas as pd
import json
from fluidfoam import readmesh
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import scienceplots

if __name__ == "__main__":

    plt.style.use('science')
    plt.rcParams.update({'font.size': 22})

    data_with = pd.read_csv("injectivity_withBHP.csv")
    data_no = pd.read_csv("injectivity_noBHP.csv")

    dx = 1          # m
    dh = 0.1        # m
    dv = dx*dh**2   # m^3
    Nx = 1000
    phi = 0.1593
    PV = Nx*dx * dh**2 * phi # m^3

    # PLOT INJECTIVITY ==================================
    
    BHP_init = 6.3e6
    gas_bkt_noBHP = 5.7e7
    gas_bkt_BHP = 4.1e7
    
    # fig, ax = plt.subplots(1,4,
    #                        figsize=(18,5),
    #                        gridspec_kw={'width_ratios': [1, 1, 1, 1]})
    fig, ax = plt.subplots(
        2, 2,
        figsize=(18, 10),
        gridspec_kw={'width_ratios': [1, 1]}
    )
    ax = ax.ravel()

    # plots with BHP limit
    ax[0].plot(data_with['t'], data_with['qt_inj']*dv, color='tab:orange', lw=2, label = 'with BHP limit')
    ax[1].plot(data_with['t'], data_with['J']*dv, color='tab:blue', lw=2, label = 'with BHP limit')
    ax[2].plot(data_with['t'], data_with['p_bh_inj'], color='tab:green', lw=2, label = 'with BHP limit')
    ax[3].plot(data_with['t'], data_with['Q_prod'], color='tab:red', lw=2, label = 'with BHP limit')

    # plots without BHP limit
    ax[0].plot(data_no['t'], data_no['qt_inj']*dv, color='tab:orange', lw=2, ls='--', label = 'without BHP limit')
    ax[1].plot(data_no['t'], data_no['J']*dv, color='tab:blue', lw=2, ls='--', label = 'without BHP limit')
    ax[2].plot(data_no['t'], data_no['p_bh_inj'], color='tab:green', lw=2, ls='--', label = 'without BHP limit')    
    ax[3].plot(data_no['t'], data_no['Q_prod'], color='tab:red', lw=2, ls='--', label = 'without BHP limit')

    ax[0].axvline(BHP_init, color='k', ls='--', lw=1.5, label = "BHP control begin")
    ax[1].axvline(BHP_init, color='k', ls='--', lw=1.5, label = "BHP control begin")
    ax[2].axvline(BHP_init, color='k', ls='--', lw=1.5, label = "BHP control begin")
    ax[0].axvline(gas_bkt_noBHP, color='gray', ls='--', lw=1.5, label = "gas breakthrough (with BHP limit)")
    ax[2].axvline(gas_bkt_BHP, color='gray', ls='--', lw=1.5, label = "gas breakthrough (no BHP limit)")
    ax[3].axvline(gas_bkt_noBHP, color='k', ls='--', lw=1.5, label = "gas breakthrough (with BHP limit)")
    ax[3].axvline(gas_bkt_BHP, color='gray', ls='--', lw=1.5, label = "gas breakthrough (no BHP limit)")

    ax[0].grid()
    ax[0].set_xlabel(r"$t$ [s]")
    ax[0].set_ylabel(r"$q_\textrm{inj}$ [$\textrm{m}^3/\textrm{s}$]")

    ax[1].grid()
    ax[1].set_xlabel(r"$t$ [s]")
    ax[1].set_ylabel(r"$J$ [$\textrm{m}^3/\textrm{s}\cdot\textrm{Pa}$]")

    ax[2].grid()
    ax[2].set_xlabel(r"$t$ [s]")
    ax[2].set_ylabel(r"$p_\textrm{bh}$ [Pa]")

    ax[3].grid()
    ax[3].set_xlabel(r"$t$ [s]")
    ax[3].set_ylabel(r"$Q_\textrm{prod}$ [$\textrm{PV}$]")

    ax[0].legend(loc='best', fontsize=14)
    ax[1].legend(loc='best', fontsize=14)
    ax[2].legend(loc='best', fontsize=14)
    ax[3].legend(loc='best', fontsize=14)

    ax[0].set_title("Injection rate")
    ax[1].set_title("Injectivity")
    ax[2].set_title("Bottom hole pressure")
    ax[3].set_title("Water production")

    plt.tight_layout()
    plt.savefig("injectivity_comparison.png", dpi=300)
    plt.show()

    

    
            