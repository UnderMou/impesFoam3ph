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
import scienceplots

if __name__ == "__main__":

    plt.style.use('science')
    plt.rcParams.update({'font.size': 22})

    # User defined ===========
    exp_name = "CTscanSA"

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

    n_folders = 1000
    # ========================

    pathCasePrep = f"casesPrep/{exp_name}"
    
    with open("exps_setup.json", "r") as f:
        exp_configs = json.load(f)

    exp_config = exp_configs[exp_name]

    # Reading the simulation results    =====================================

    # mesh
    x, y, z = readmesh('templates/CTscan', structured=True)
    z = z.flatten()
    L = z[-1]
    xD = z/L
    # print(z)

    # time
    ppt = exp_config['postProcTime']
    start, end, step = ppt['start'], ppt['end'], ppt['step']

    t = np.arange(start, end, step)
    print(t)
    U = 1.61572e-5
    phi = 0.22
    Swc = 0.197
    Sgr = 0.013
    Sor = 0.103
    tD = U * t / (L*phi*(1-Swc-Sgr-Sor))
    print(tD)
    

    t_des = 1000
    tdes_id = np.argmin(np.abs(t - t_des))
    print(t[tdes_id])

    # # results
    # n_folders = sum(
    #     os.path.isdir(os.path.join(exp_config['output_path'], item))
    #     for item in os.listdir(exp_config['output_path'])
    # )
    # print(n_folders)

    # PLOT SO, MRF, and Fstars  =========================
    if False:

        fig, ax = plt.subplots(5, 1, figsize=(8, 12))
        for i in range(n_folders):
            exp_dir = exp_config['output_path'] + "/run_" + str(i)
            data_single = uqt.parse_openfoam_case(
                case_dir=exp_dir,
                variables=["Sa", "Sb", "Sc", "MRF", "Fdry", "Fshear", "Foil"],
                time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
            )

            Sc = data_single.Sc[tdes_id,:]
            MRF = data_single.MRF[tdes_id,:]
            Fdry = data_single.Fdry[tdes_id,:]
            Fshear = data_single.Fshear[tdes_id,:]
            Foil = data_single.Foil[tdes_id,:]

            ax[0].plot(xD, Sc, c='gray',alpha=0.4, linewidth=1.0)
            ax[1].plot(xD, 1/MRF, c='gray',alpha=0.4, linewidth=1.0)
            ax[2].plot(xD, Fdry, c='gray',alpha=0.4, linewidth=1.0)
            ax[3].plot(xD, Fshear, c='gray',alpha=0.4, linewidth=1.0)
            ax[4].plot(xD, Foil, c='gray',alpha=0.4, linewidth=1.0)

            print(i,f" t={data_single['time'][tdes_id].item()} | ", n_folders)

        ax[1].set_ylim([0,1000])
        ax[1].axhline(1000.0, c='k', ls='--', label='MRF=1000')
        
        ax[0].set_ylabel(r"$S_o$ [-]")
        ax[1].set_ylabel(r"MRF [-]")
        ax[2].set_ylabel(r"$F_{dry}$ [-]")
        ax[3].set_ylabel(r"$F_{shear}$ [-]")
        ax[4].set_ylabel(r"$F_{oil}$ [-]")
        ax[4].set_xlabel(r"$x_D$ [-]")

        plt.grid()
        ax[1].legend(loc='best')
        plt.savefig(exp_config['output_path']+"/stars_profiles.pdf", dpi=300)
        plt.show()
        plt.close()

    # PLOT SO and Samples name  =========================
    
    if True:

        fig, ax = plt.subplots(1, 1, figsize=(7, 5))
        for i in range(n_folders):
            exp_dir = exp_config['output_path'] + "/run_" + str(i)
            data_single = uqt.parse_openfoam_case(
                case_dir=exp_dir,
                variables=["Sc"],
                time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
            )

            Sc = data_single.Sc[tdes_id,:]

            ax.plot(xD, Sc, c='gray',alpha=0.4, linewidth=1.5)

            print(i,f" t={data_single['time'][tdes_id].item()} | ", n_folders)

        for sample_name, idx in dictSamples.items():
            exp_dir = exp_config['output_path'] + "/run_" + str(idx)
            data_single = uqt.parse_openfoam_case(
                case_dir=exp_dir,
                variables=["Sc"],
                time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
            )
            Sc = data_single.Sc[tdes_id,:]

            color = sample_colors[sample_name]

            plt.plot(xD, Sc,color=color,linewidth=3.0,label=f'Sample {sample_name}')
        
        ax.set_ylabel(r"$S_o$ [-]")
        ax.set_xlabel(r"$x_D$ [-]")

        plt.grid()
        ax.legend(loc='best')
        plt.tight_layout()
        plt.savefig(exp_config['output_path']+"/samples_profiles.pdf", dpi=300)
        plt.show()
        plt.close()


    # Plot phases saturation for all samples at desired time =========================

    if True:

        t_des = [1100, 2450] # [1000, 2350]  
        tD_des = [0.44, 1.0]

        fig, ax = plt.subplots(2, 3, figsize=(16, 9), sharex=True)

        for j in range(len(t_des)):

            tdes_id = np.argmin(np.abs(t - t_des[j]))

            for i in range(n_folders):
                exp_dir = exp_config['output_path'] + "/run_" + str(i)
                data_single = uqt.parse_openfoam_case(
                    case_dir=exp_dir,
                    variables=["Sa", "Sb", "Sc"],
                    time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
                )

                Sc = data_single.Sc[tdes_id,:]
                Sb = data_single.Sb[tdes_id,:]
                Sa = data_single.Sa[tdes_id,:]

                ax[j,0].plot(xD, Sc, c='gray',alpha=0.4, linewidth=1.0)
                ax[j,1].plot(xD, Sb, c='gray',alpha=0.4, linewidth=1.0)
                ax[j,2].plot(xD, Sa, c='gray',alpha=0.4, linewidth=1.0)

                if j == 0 and i == 0:
                    ax[0,0].plot(xD, Sc, c='gray',alpha=0.4, linewidth=1.0, label='Samples')
                
                ax[j,1].set_title(f"$t_D={tD_des[j]} [PVI]$")

                print(i,f" t={data_single['time'][tdes_id].item()} | ", n_folders)

            # experimental plot
            df = pd.read_csv(f"{exp_config['input_path']}/experimental_results/So_PVI_{tD_des[j]}.csv")
            ax[j,0].plot(df['x'], df[' y'], color='k', linestyle='-', marker='o', markersize=5, label='CT scan (Tang, J., 2019)')
            df = pd.read_csv(f"{exp_config['input_path']}/experimental_results/Sw_PVI_{tD_des[j]}.csv")
            ax[j,1].plot(df['x'], df[' y'], color='k', linestyle='-', marker='o', markersize=5)
            df = pd.read_csv(f"{exp_config['input_path']}/experimental_results/Sg_PVI_{tD_des[j]}.csv")
            ax[j,2].plot(df['x'], df[' y'], color='k', linestyle='-', marker='o', markersize=5)

        ax[0,0].set_ylabel(r"$S_o$ [-]")
        ax[0,1].set_ylabel(r"$S_w$ [-]")
        ax[0,2].set_ylabel(r"$S_g$ [-]")
        ax[1,0].set_ylabel(r"$S_o$ [-]")
        ax[1,1].set_ylabel(r"$S_w$ [-]")
        ax[1,2].set_ylabel(r"$S_g$ [-]")
        
        ax[1,0].set_xlabel(r"$x/L$ [-]")
        ax[1,1].set_xlabel(r"$x/L$ [-]")
        ax[1,2].set_xlabel(r"$x/L$ [-]")
        

        ax[0,0].legend(loc='best')

        for a in ax.flat:
            a.grid(True, alpha=0.5)
            a.set_ylim([0, 1])
        
        fig.subplots_adjust(
            wspace=0.3,
            hspace=0.40
        )

        ax[0,0].legend(
            frameon=False,
            loc='upper center',
            bbox_to_anchor=(2.7, 1.4),
            ncol=2
        )

        ax[0,0].text(
            -0.2, 1.2, '(a)',
            transform=ax[0,0].transAxes,
            ha='left', va='top',
            fontsize=22,
            fontweight='bold'
        )

        ax[1,0].text(
            -0.2, 1.2, '(b)',
            transform=ax[1,0].transAxes,
            ha='left', va='top',
            fontsize=22,
            fontweight='bold'
        )

        plt.savefig(exp_config['output_path']+"/samples_So.pdf", dpi=300)
        plt.show()
        plt.close()
