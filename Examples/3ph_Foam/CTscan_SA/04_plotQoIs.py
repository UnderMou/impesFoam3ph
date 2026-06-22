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

    # User defined ===========
    exp_name = "CTscanSA"

    # CTscan info
    A = 0.001256637 # [m2] 
    L = 0.4         # [m]  
    phi = 0.22      # [-]  

    # Samples dict
    dictSamples = {
        'S1': 7, 
        'S2': 4,
        'S3': 13
    }
    sample_colors = {
        'S1': 'r',
        'S2': 'b',
        'S3': 'k',
    }

    # ========================

    with open("exps_setup.json", "r") as f:
        exp_configs = json.load(f)

    exp_config = exp_configs[exp_name]

    df = pd.read_csv('casesPrep/CTscanSA/input_samples.csv')

    df = df.iloc[:1000]

    ###########################

    ppt = exp_config['postProcTime']
    # start, end, step = ppt['start'], ppt['end'], ppt['step']
    start, end, step = ppt['start'], 10000, ppt['step'] # TODO

    # PLOT CUMULATIVE PHASE PRODUCTION FOR ALL SAMPLES ==================================

    qo = np.zeros((df.shape[0], int((end-start)/(step))))
    qw = np.zeros((df.shape[0], int((end-start)/(step))))
    qg = np.zeros((df.shape[0], int((end-start)/(step))))

    for ind,i in enumerate(df.index):
        print(i,ind)

        exp_dir = exp_config['output_path'] + "/run_" + str(i)
        data_single = uqt.parse_openfoam_case(
            case_dir=exp_dir,
            variables=["Sa", "Sb", "Sc", "Fa", "Fb", "U"],
            time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
        )

        # Check along all simulation time
        for j in range(len(data_single.time)):
            U = np.asarray(data_single.U[j,:])
            Fa = np.asarray(data_single.Fa[j,:])
            Fb = np.asarray(data_single.Fb[j,:])
            Fc = 1.0 - Fa - Fb
            Sc = np.asarray(data_single.Sc[j,:])

            outid = -2
            Uo_out = U[outid][2] * Fc[outid]    # U is the z component
            Uw_out = U[outid][2] * Fb[outid]    # U is the z component
            Ug_out = U[outid][2] * Fa[outid]    # U is the z component
            qoil = Uo_out * A             # m3/s
            qwater = Uw_out * A           # m3/s
            qgas = Ug_out * A             # m3/s

            qo[ind,j] = qoil
            qw[ind,j] = qwater
            qg[ind,j] = qgas
   
    # eval through time
    t = np.asarray(data_single.time)
    dt = t[1:] - t[:-1]

    Qoil = np.zeros((df.shape[0], int((end-start)/(step))-1))
    Qwater = np.zeros((df.shape[0], int((end-start)/(step))-1))
    Qgas = np.zeros((df.shape[0], int((end-start)/(step))-1))

    PV = L*A*phi
    OOIP = 0.4599

    Uin = U[-1][2]
    PVI = Uin * t / (phi*L*(1 - 0.197 - 0.103 - 0.013))
    
    # Cumulative phase production
    fig, ax = plt.subplots(1,2,
                           figsize=(12,5),
                           gridspec_kw={'width_ratios': [1, 1.25]})

    # Normalize swR_s values to the colormap range
    norm = mcolors.Normalize(vmin=min(df.sfdry), vmax=max(df.sfdry))
    cmap = cm.viridis  # choose any colormap you like

    for i in range(Qoil.shape[0]):
        for j in range(Qoil.shape[1]):
            Qoil[i,j] = np.sum(qo[i,:j]*dt[:j])
            Qwater[i,j] = np.sum(qw[i,:j]*dt[:j])
            Qgas[i,j] = np.sum(qg[i,:j]*dt[:j])

        color = cmap(norm(df.sfdry.iloc[i]))

        ax[0].plot(PVI[1:],(Qoil[i,:]/PV),color=color, linewidth=1.5)
        ax[1].plot(PVI[1:],(Qgas[i,:]/PV),color=color, linewidth=1.5)

    for sample_name, idx in dictSamples.items():
        color = sample_colors[sample_name]
        ax[0].plot(PVI[1:],Qoil[idx, :] / PV,color=color,linewidth=2.5,label=f'Sample {sample_name}')
        ax[1].plot(PVI[1:],Qgas[idx, :] / PV,color=color,linewidth=2.5,label=f'Sample {sample_name}')
    
    # Create colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax[1])
    cbar.set_label(r'$SF$')

    ax[0].grid()
    ax[0].set_xlabel(r"$t_D$ [PVI]")
    ax[0].set_ylabel(r"$Q_o$ [PV]")

    ax[1].grid()
    ax[1].set_xlabel(r"$t_D$ [PVI]")
    ax[1].set_ylabel(r"$Q_g$ [PV]")

    ax[0].text(
        -0.2, 1.2, '(a)',
        transform=ax[0].transAxes,
        ha='left', va='top',
        fontsize=22,
        fontweight='bold'
    )

    ax[1].text(
        -0.2, 1.2, '(b)',
        transform=ax[1].transAxes,
        ha='left', va='top',
        fontsize=22,
        fontweight='bold'
    )
    
    
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(exp_config['output_path'] + "/cumulativeProd.pdf", dpi=300)
    plt.show()

    # PLOT ORF AND GOR FOR ALL SAMPLES ==================================

    fig, ax = plt.subplots(1,2,
                           figsize=(12,5),
                           gridspec_kw={'width_ratios': [1, 1.25]})

    # Normalize swR_s values to the colormap range
    norm = mcolors.Normalize(vmin=min(df.sfdry), vmax=max(df.sfdry))
    cmap = cm.viridis  # choose any colormap you like

    for i in range(Qoil.shape[0]):

        color = cmap(norm(df.sfdry.iloc[i]))

        ax[0].plot(PVI[1:],((Qoil[i,:]/PV)/OOIP)*100,color=color, linewidth=1.5)
        ax[1].plot(PVI[1:],qg[i,1:]/qo[i,1:],color=color, linewidth=1.5)

    for sample_name, idx in dictSamples.items():
        color = sample_colors[sample_name]
        ax[0].plot(PVI[1:],((Qoil[idx,:]/PV)/OOIP)*100,color=color,linewidth=2.5,label=f'Sample {sample_name}')
        ax[1].plot(PVI[1:],qg[idx,1:]/qo[idx,1:],color=color,linewidth=2.5,label=f'Sample {sample_name}')
    
    # Create colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax[1])
    cbar.set_label(r'$SF$')

    ax[0].grid()
    ax[0].set_xlabel(r"$t_D$ [PVI]")
    ax[0].set_ylabel(r"ORF [\%]")
    ax[0].set_ylim([0,60])

    ax[1].grid()
    ax[1].set_xlabel(r"$t_D$ [PVI]")
    ax[1].set_ylabel(r"GOR [$q_g/q_o$]")
    
    ax[0].text(
        -0.2, 1.2, '(c)',
        transform=ax[0].transAxes,
        ha='left', va='top',
        fontsize=22,
        fontweight='bold'
    )

    ax[1].text(
        -0.2, 1.2, '(d)',
        transform=ax[1].transAxes,
        ha='left', va='top',
        fontsize=22,
        fontweight='bold'
    )

    
    # plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(exp_config['output_path'] + "/orfgor.pdf", dpi=300)
    plt.show()       

            