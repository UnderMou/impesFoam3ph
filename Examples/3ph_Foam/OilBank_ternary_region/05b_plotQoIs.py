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
    exp_name = "CTscan"

    # CTscan info
    A = 0.001256637 # [m2] 
    L = 0.4         # [m]  
    phi = 0.22      # [-]  

    # ========================

    with open("exps_setup.json", "r") as f:
        exp_configs = json.load(f)

    exp_config = exp_configs[exp_name]

    Xfile = 'X.csv'

    df = pd.read_csv(exp_config['output_path'] + '/' + Xfile)
    filename = Xfile.split('.')[0]

    ###########################

    # Get R and L states as criteria bellow:
    
    # tolerance
    tol = 1e-2

    # condition
    # SoR = 0.4 and SgR = 0.013
    so = 1 - df["sw"] - df["sg"]
    mask = (
        np.isclose(df["sg"], 0.03, atol=tol) &
        np.isclose(so, 0.4, atol=tol)
    )

    # indexes
    indexes = df.index[mask]
    # indexes = indexes[:3]
    
    print(indexes)
    print(df.iloc[indexes])

    ppt = exp_config['postProcTime']
    # start, end, step = ppt['start'], ppt['end'], ppt['step']
    start, end, step = ppt['start'], 50000, ppt['step'] # TODO

    # plot oil saturation profiles at desired time
    t_des = 15000 # 9000
    tempMesh = exp_config['input_path']

    x, y, z = readmesh(tempMesh, structured=True)
    z = np.asarray(z).flatten()
    nx = z.shape[0]
    L = z[-1]
    
    Scs = np.zeros((indexes.shape[0], nx))
    swL_s = np.zeros((indexes.shape[0], 1))

    for ind,i in enumerate(indexes):
        print(i,ind)

        exp_dir = exp_config['output_path'] + "/run_" + str(i)
        data_single = uqt.parse_openfoam_case(
            case_dir=exp_dir,
            variables=["Sc"],
            time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
        )

        # Check along all simulation time
        for j in range(len(data_single.time)):

            if data_single.time[j] == t_des:

                Scs[ind,:] = np.asarray(data_single.Sc[j,:])
                swL_s[ind] = df["sw_inj"][indexes[ind]]
    

    # PLOT OIL SATURATION PROFILES  =========================================

    fmoil = 0.823
    floil = 0.295
    xD = z/L

    fig, ax = plt.subplots(1,1,
                           figsize=(7,7),
                           gridspec_kw={'width_ratios': [1.25]})

    # Normalize swR_s values to the colormap range
    norm = mcolors.Normalize(vmin=min(swL_s), vmax=max(swL_s))
    cmap = cm.viridis  # choose any colormap you like


    for i in range(len(indexes)):
        color = cmap(norm(swL_s[i]))
        ax.plot(xD, Scs[i, :], color=color, linewidth=3)
    
    ax.axhline(fmoil, color='k', linestyle='--', linewidth=1.0)
    ax.axhline(floil, color='k', linestyle='--', linewidth=1.0)
    x_text = xD.min()
    ax.text(x_text, fmoil + 0.001, 'fmoil',
            ha='left', va='bottom', fontsize=16)

    ax.text(x_text, floil + 0.001, 'floil',
            ha='left', va='bottom', fontsize=16)
    
    # Create colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r'$S_{w}^R$')

    ax.set_xlabel(r"$x/L$")
    ax.set_ylabel(r"$S_o$")
    ax.grid()

    plt.savefig(exp_config['output_path'] + "/oilProfiles.pdf", dpi=300)
    plt.show()
    plt.close()



    # PLOT CUMULATIVE PHASE PRODUCTION FOR ALL SAMPLES ==================================

    qo = np.zeros((indexes.shape[0], int((end-start)/(step))))
    qw = np.zeros((indexes.shape[0], int((end-start)/(step))))
    qg = np.zeros((indexes.shape[0], int((end-start)/(step))))

    for ind,i in enumerate(indexes):
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
          
            Uo_out = U[-1][2] * Fc[-1]    # U is the z component
            Uw_out = U[-1][2] * Fb[-1]    # U is the z component
            Ug_out = U[-1][2] * Fa[-1]    # U is the z component
            qoil = Uo_out * A             # m3/s
            qwater = Uw_out * A           # m3/s
            qgas = Ug_out * A             # m3/s

            qo[ind,j] = qoil
            qw[ind,j] = qwater
            qg[ind,j] = qgas
   
    # eval through time
    t = np.asarray(data_single.time)
    dt = t[1:] - t[:-1]

    Qoil = np.zeros((indexes.shape[0], int((end-start)/(step))-1))
    Qwater = np.zeros((indexes.shape[0], int((end-start)/(step))-1))
    Qgas = np.zeros((indexes.shape[0], int((end-start)/(step))-1))

    PV = L*A*phi
    OOIP = 0.4

    Uin = U[-1][2]
    PVI = Uin * t / (phi*L*(1 - 0.197 - 0.103 - 0.013))
    
    # Cumulative phase production
    fig, ax = plt.subplots(1,2,
                           figsize=(12,5),
                           gridspec_kw={'width_ratios': [1, 1.25]})

    # Normalize swR_s values to the colormap range
    norm = mcolors.Normalize(vmin=min(swL_s), vmax=max(swL_s))
    cmap = cm.viridis  # choose any colormap you like

    for i in range(Qoil.shape[0]):
        for j in range(Qoil.shape[1]):
            Qoil[i,j] = np.sum(qo[i,:j]*dt[:j])
            Qwater[i,j] = np.sum(qw[i,:j]*dt[:j])
            Qgas[i,j] = np.sum(qg[i,:j]*dt[:j])

        color = cmap(norm(swL_s[i]))

        ax[0].plot(PVI[1:],(Qoil[i,:]/PV),color=color, linewidth=3)
        ax[1].plot(PVI[1:],(Qgas[i,:]/PV),color=color, linewidth=3)
    
    # Create colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax[1])
    cbar.set_label(r'$S_{w}^R$')

    ax[0].grid()
    ax[0].set_xlabel(r"$t_D$ [PVI]")
    ax[0].set_ylabel(r"$Q_o$ [PV]")

    ax[1].grid()
    ax[1].set_xlabel(r"$t_D$ [PVI]")
    ax[1].set_ylabel(r"$Q_g$ [PV]")
    
    plt.savefig(exp_config['output_path'] + "/cumulativeProd.pdf", dpi=300)
    plt.tight_layout()
    plt.show()

    # PLOT Oil-cut AND Gas-cut FOR ALL SAMPLES ==================================

    fig, ax = plt.subplots(1,2,
                           figsize=(12,5),
                           gridspec_kw={'width_ratios': [1, 1.25]})

    # Normalize swR_s values to the colormap range
    norm = mcolors.Normalize(vmin=min(swL_s), vmax=max(swL_s))
    cmap = cm.viridis  # choose any colormap you like

    for i in range(Qoil.shape[0]):

        color = cmap(norm(swL_s[i]))

        ax[0].plot(t, qo[i,:]/PV,color=color, linewidth=3)
        ax[1].plot(t, qg[i,:]/PV,color=color, linewidth=3)
    
    # Create colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax[1])
    cbar.set_label(r'$S_{w}^R$')

    ax[0].grid()
    ax[0].set_xlabel(r"$t$ [sec]")
    ax[0].set_ylabel(r"Oil-cut [PV/sec]")

    ax[1].grid()
    ax[1].set_ylim([0,3e-5])
    ax[1].set_xlabel(r"$t$ [sec]")
    ax[1].set_ylabel(r"Gas-cut [PV/sec]")
    ax[1].ticklabel_format(axis='y',
                    style='sci',
                    scilimits=(0, 0),
                    useMathText=True)
    

    plt.savefig(exp_config['output_path'] + "/oilGas_cut.pdf", dpi=300)
    plt.tight_layout()
    plt.show()    

    # PLOT ORF AND GOR FOR ALL SAMPLES ==================================

    fig, ax = plt.subplots(1,2,
                           figsize=(12,5),
                           gridspec_kw={'width_ratios': [1, 1.25]})

    # Normalize swR_s values to the colormap range
    norm = mcolors.Normalize(vmin=min(swL_s), vmax=max(swL_s))
    cmap = cm.viridis  # choose any colormap you like

    for i in range(Qoil.shape[0]):

        color = cmap(norm(swL_s[i]))

        ax[0].plot(PVI[1:],((Qoil[i,:]/PV)/OOIP)*100,color=color, linewidth=3)
        ax[1].plot(PVI[1:],qg[i,1:]/qo[i,1:],color=color, linewidth=3)
    
    # Create colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax[1])
    cbar.set_label(r'$S_{w}^R$')

    ax[0].grid()
    ax[0].set_xlabel(r"$t_D$ [PVI]")
    ax[0].set_ylabel(r"ORF [\%]")

    ax[1].grid()
    ax[1].set_xlabel(r"$t_D$ [PVI]")
    ax[1].set_ylabel(r"GOR [$q_g/q_o$]")
    

    plt.savefig(exp_config['output_path'] + "/orfgor.pdf", dpi=300)
    plt.tight_layout()
    plt.show()       

            