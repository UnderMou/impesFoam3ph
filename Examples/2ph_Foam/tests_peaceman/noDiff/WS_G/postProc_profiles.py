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

    exp_dir = "./"
    data_single = uqt.parse_openfoam_case(
        case_dir=exp_dir,
        variables=["qt","qb","p","p_bh","WI","mob_t","Cs","Sb"],
        time_dirs=[f"{i:g}" for i in np.arange(1e7, 5e7, 1e6)]   
        # time_dirs=[f"{i:g}" for i in np.arange(1e6, 5e7, 1e6)]   
    )    

    # PLOT INJECTIVITY ==================================

    t = np.asarray(data_single.time)
    
    bkt = 5.7e7
    max_BHP = 2.0e7
    saveName = 'injectivity_withBHP'

    # bkt = 4.1e7
    # max_BHP = 1.0e8
    # saveName = 'injectivity_noBHP'
    
    dx = 1          # m
    dh = 0.1        # m
    dv = dx*dh**2   # m^3
    Nx = 1000
    phi = 0.1593
    PV = Nx*dx * dh**2 * phi # m^3

    J = np.zeros_like(t)
    qt_inj = np.zeros_like(t)
    qw_prod = np.zeros_like(t)
    p_bh_inj = np.zeros_like(t)
    p_inj = np.zeros_like(t)
    area_Csw = np.zeros_like(t)
    for j in range(len(data_single.time)):
        qt = np.asarray(data_single.qt[j,:])
        qb = np.asarray(data_single.qb[j,:])
        p_bh = np.asarray(data_single.p_bh[j,:])
        WI = np.asarray(data_single.WI[j,:])
        mob_t = np.asarray(data_single.mob_t[j,:])
        p = np.asarray(data_single.p[j,:])
        Cs = np.asarray(data_single.Cs[j,:])
        Sb = np.asarray(data_single.Sb[j,:])

        # domain arrays
        qt_inj[j] = qt[0]
        p_bh_inj[j] = p_bh[0]
        p_inj[j] = p[0]
        qw_prod[j] = -qb[-1]*dv

        J[j] = WI[0]*mob_t[0]
        # J[j] = qt_inj[j]/(p_bh_inj[j] - p_inj[j])

        area_Csw[j] = np.sum(dx*Cs*Sb)

    Q_prod = np.zeros(len(qw_prod))
    Q_prod[0] = 0
    dt = t[1:] - t[:-1]
    for j in range(Q_prod[1:].shape[0]+1):
        Q_prod[j] = np.sum(qw_prod[:j]*dt[:j])
    
    fig, ax = plt.subplots(1,4,
                           figsize=(18,5),
                           gridspec_kw={'width_ratios': [1, 1, 1, 1]})

    ax[0].plot(t, qt_inj, color='tab:orange', lw=2)
    ax[1].plot(t, J, color='tab:blue', lw=2, label=r"$J=\frac{q_t}{p_\textrm{bh}-p_\textrm{res. cell}}$")
    ax[2].plot(t, p_bh_inj, color='tab:green', lw=2, label = "BHP")
    # ax[2].plot(t, p_inj, color='tab:purple', lw=2, label = r"$p_{\textrm{res. cell}}$")
    ax[3].plot(t, Q_prod/PV, color='tab:red', lw=2, label = "prod rate")

    # ax[0].axvline(bkt, color='k', ls='--', lw=1.5, label = "gas breakthrough")
    # ax[1].axvline(bkt, color='k', ls='--', lw=1.5)
    # ax[2].axvline(bkt, color='k', ls='--', lw=1.5)
    # ax[2].axhline(max_BHP, color='r', ls='--', lw=1, label="max BHP")
    # ax[3].axvline(bkt, color='k', ls='--', lw=1.5)

    ax[0].grid()
    ax[0].set_xlabel(r"$t$ [s]")
    ax[0].set_ylabel(r"$q_\textrm{inj}$ [$\textrm{m}^3/\textrm{s}$]")
    ax[0].set_ylim([0, 1.1e-6])

    ax[1].grid()
    ax[1].set_xlabel(r"$t$ [s]")
    ax[1].set_ylabel(r"$J$ [$\textrm{m}^3/\textrm{s}\cdot\textrm{Pa}$]")
    ax[1].set_ylim([0, 1.5e-9])

    ax[2].grid()
    ax[2].set_xlabel(r"$t$ [s]")
    ax[2].set_ylabel(r"$p_\textrm{bh}$ [Pa]")

    ax[3].grid()
    ax[3].set_xlabel(r"$t$ [s]")
    ax[3].set_ylabel(r"$Q_\textrm{prod}$ [$\textrm{PV}$]")
    ax[3].set_ylim([-1e-2, 0.1])

    ax[0].legend(loc='best', fontsize=14)
    ax[1].legend(loc='best', fontsize=24)
    ax[2].legend(loc='best', fontsize=14)

    ax[0].set_title("Injection rate")
    ax[1].set_title("Injectivity")
    ax[2].set_title("Bottom hole pressure")
    ax[3].set_title("Water production")

    plt.tight_layout()
    plt.savefig(saveName + ".png", dpi=300)
    plt.show()
    plt.close()

    save_injec = False
    if save_injec:
        df = pd.DataFrame({'t': t, 'qt_inj': qt_inj, 'J': J, 'p_bh_inj': p_bh_inj, 'Q_prod': Q_prod/PV})
        df.to_csv(saveName + '.csv', index=False)

    
    plt.figure(figsize=(6,5))
    plt.plot(t,p_bh_inj-p_inj)
    plt.xlabel(r"$t$ [s]")
    plt.ylabel(r"$p_\textrm{bh}-p_\textrm{res. cell}$ [Pa]")
    plt.title(r"$\Delta p =$ BHP - reservoir pressure")
    plt.grid()
    plt.tight_layout()
    plt.savefig(saveName + "_deltaP.png", dpi=300)
    plt.show()

    plt.figure(figsize=(6,5))
    plt.plot(t,area_Csw)
    plt.xlabel(r"$t$ [s]")
    plt.ylabel(r"$\int_0^L C_s \cdot S_w \, dx$")
    plt.title("Surfactant area")
    plt.grid()
    plt.tight_layout()
    # plt.savefig(saveName + "_areaCsSw.png", dpi=300)
    plt.show()