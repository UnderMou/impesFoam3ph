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
    exp_name = ("Lyu_fmmobL", "Lyu", "Lyu_fmmobH")
    fmmob_values = {
        "Lyu_fmmobL": 875,
        "Lyu": 1750,
        "Lyu_fmmobH": 3500
    }

    # desired case
    SwR = 0.6
    SoR = 0.24
    SwL = 0.3
    tdes = 50000

    thres_min = -0.01
    thres_max = 0.01
    # ========================

    
    fig, ax = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    for i in range(len(exp_name)):
    
        with open("exps_setup.json", "r") as f:
            exp_configs = json.load(f)

        exp_config = exp_configs[exp_name[i]]

        # Reading the simulation results  

        n_folders = sum(
            os.path.isdir(os.path.join(exp_config['output_path'], item))
            for item in os.listdir(exp_config['output_path'])
        )

        with open(exp_config['output_path']+"/samples.csv", "r") as f:
            reader = csv.reader(f)
            next(reader) 
            data = [list(map(float, row)) for row in reader]
        X = np.array(data)
        sw = X[:,0]
        sg = X[:,1]
        so = 1 - sw - sg
        swL = X[:,2]

        id = np.argmin(np.abs(sw - SwR) + np.abs(so - SoR) + np.abs(swL - SwL))
        print(f"Case {exp_name[i]}, run_{id}: Sw={sw[id]:.3f}, Sg={sg[id]:.3f}, So={so[id]:.3f}")
  
        ppt = exp_config['postProcTime']
        start, end, step = ppt['start'], ppt['end'], ppt['step']
    
        exp_dir = exp_config['output_path'] + "/run_" + str(id)
        data_single = uqt.parse_openfoam_case(
            case_dir=exp_dir,
            variables=["Sa", "Sb", "Sc", "Fa"],
            time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
        )
        idt = np.argmin(np.abs(data_single.time.to_numpy() - tdes))

        # mesh
        x, y, z = readmesh(exp_dir, structured=True)
        x = x.flatten()
        L = x[-1]
        xD = x/L
        
        print(f"Case {exp_name[i]}, run_{id}: t={data_single.time[idt]:.3f} s")
        
        Sc = data_single.Sc[idt,:]

        ax[0].plot(xD, Sc, linewidth=2, label=f"fmmob = {fmmob_values[exp_name[i]]}")
        ax[1].plot(xD, np.gradient(Sc), linewidth=2)

    ax[1].axhline(thres_max, color='k', linewidth = 2, linestyle='--', label=r"$\partial S_o / \partial x$ thresholds")
    ax[1].axhline(thres_min, color='k',linewidth = 2,  linestyle='--')
    
    ax[1].set_xlabel(r"$x/L$ [-]")
    ax[1].set_ylabel(r"$\partial S_o / \partial x$ [-]")
    ax[0].set_ylabel(r"$S_o$ [-]")
    ax[0].legend(loc='best', fontsize=16)
    ax[1].legend(loc='best', fontsize=16)
    for a in ax:
        a.grid(True, alpha=0.75)

    ax[0].text(
        -0.1, 1.1, '(a)',
        transform=ax[0].transAxes,
        ha='left', va='top',
        fontsize=22,
        fontweight='bold'
    )

    ax[1].text(
        -0.1, 1.1, '(b)',
        transform=ax[1].transAxes,
        ha='left', va='top',
        fontsize=22,
        fontweight='bold'
    )
    
    plt.savefig(exp_config['output_path']+"/criteria.pdf", dpi=300)
    plt.show()
                


    
