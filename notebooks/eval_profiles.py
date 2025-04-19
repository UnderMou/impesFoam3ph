import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from jinja2 import Environment, FileSystemLoader
from functools import partial
from tqdm import tqdm
from multiprocessing import Pool
import itertools
import scienceplots
from matplotlib.colors import Normalize
# from matplotlib.cm import get_cmap
import matplotlib.cm as cm

from aux_funcs import *

def Fdry(Sw, params):
    _,SF,sfbet,_,_,_,_ = params
    print(SF)
    Fdry = np.zeros_like(Sw)
    for i in range(len(Sw)):
        Fdry = 0.5 + (1.0/np.pi)*np.arctan(sfbet*(Sw-SF))
    return Fdry

def tFind(tD):
    if tD in [.36, .44]:  
        t_find = 800 #750
        isCTscan = True
        isDARTS = True
    elif tD in [.82, 1.0]:  
        t_find = 2600#2200
        isCTscan = True
        isDARTS = True
    elif tD in [.6]:
        t_find = 1400 
        isCTscan = False
        isDARTS = False
    elif tD == 10.40: 
        raise KeyError("The 'tD' provided is not expected")
        t_find = 1
        isCTscan = True
        isDARTS = False
    else: #raise KeyError("The 'tD' provided is not expected")
        raise KeyError("The 'tD' provided is not expected")
        t_find = 10000 # 39800
        isCTscan = False
        isDARTS = False
    return t_find, isCTscan, isDARTS

def checkOscilations(S):
    tol = 1e-4
    # print(S[1:])
    # print(S[:-1])
    boo1 = S[1:] < S[:-1]
    # boo2 = np.abs(S[1:] - S[:-1]) < tol
    # print(boo1,'\n',boo2)
    boo2 = boo1
    count = 0
    last = boo2[0]
    for i in range(len(boo2[1:])):
        if boo2[i] != last:
            count += 1
            last = boo2[i]
        
    
    return count

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams.update({'font.size': 16})

    tD = 0.82
    # experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_fitted_v2'
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_fitted_fixedSF'

    # Residual saturations:
    residuals = {
        'Swc': 0, # 0.197,
        'Sgr': 0, # 0.013,
        'Sor': 0  # 0.103
    }

    # Time controls:
    ti = 0
    tf = 10000 
    write_interval = 200
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    t_find, isCTscan, isDARTS = tFind(tD)
    idx = np.argmin(np.abs(t-t_find))
    

    # CT scan experimental files:
    if isCTscan:
        Sg_CT_file = 'Tang/Sg_PVI_'+str(tD)+'.csv'
        Sw_CT_file = 'Tang/Sw_PVI_'+str(tD)+'.csv'
        So_CT_file = 'Tang/So_PVI_'+str(tD)+'.csv'
        CTscan_results = get_CTscanData(So_CT_file, Sg_CT_file, Sw_CT_file)

    # DARTS files:
    if isDARTS:
        Sg_DARTS_file = 'Tang/Sg_PVI_'+str(tD)+'_DARTS.csv'
        Sw_DARTS_file = 'Tang/Sw_PVI_'+str(tD)+'_DARTS.csv'
        So_DARTS_file = 'Tang/So_PVI_'+str(tD)+'_DARTS.csv'
        DARTS_results = get_DARTSdata(Sg_DARTS_file, Sw_DARTS_file, So_DARTS_file)

    # Sampling results:
    experiment_dir = set_experimentDir(experiment_name)
    folders = get_folders(experiment_dir)
    print(folders)
    # folders = folders[:470]
    # indexes_to_remove = [55, 63, 65, 66, 136, 137, 175, 184, 205, 226, 232, 247, 311, 381, 384, 433]    # original
    folders = folders[:500]
    indexes_to_remove = [18, 19, 21, 36, 49, 57, 66, 182, 222, 236, 245, 284, 305, 310, 340, 388, 402, 412, 416, 456, 482, 493, 494] # fixedSF
    folders = remove_indexes(folders, indexes_to_remove) 
    Sas, Sbs, Scs = get_Saturations(folders, experiment_dir, idx, residuals)
    Fshears, Foils = get_FsSTARS(folders, experiment_dir, idx, residuals)
    print("Fshears:", Fshears)
    print("Type of Fshears:", type(Fshears))
    # exit()
    X_ED = pd.DataFrame(get_samples(experiment_dir))
    X_ED = X_ED.drop(indexes_to_remove).reset_index(drop=True)
    nSamples = len(folders)
    X_ED = X_ED[:nSamples]
    nCells = len(Sas[0])
    x = np.linspace(0,1,nCells)
    L = x[-1]
    x = x/L

    # Check oscilations
    nAvoid = -2
    indexes_to_remove = []
    for i,sample in enumerate(folders):
        count = checkOscilations(Scs[i][:nAvoid])
        if count >= 5:
            indexes_to_remove.append(i)
        print(i, count)

    print(indexes_to_remove)
    # folders = remove_indexes(folders, indexes_to_remove) 
    # Sas = remove_indexes(Sas, indexes_to_remove) 
    # Sbs = remove_indexes(Sbs, indexes_to_remove) 
    # Scs = remove_indexes(Scs, indexes_to_remove) 

    np.savetxt(experiment_dir + '/Scs_' + str(tD) + '.csv', Scs)

    # Plotting:
    # nAvoid = -1
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    # Samples profile
    axes[0].plot(x[:nAvoid],Scs[0][:nAvoid],c='gray',alpha=0.3,label='Samples')
    axes[1].plot(x[:nAvoid],Sbs[0][:nAvoid],c='gray',alpha=0.3,label='Samples')
    axes[2].plot(x[:nAvoid],Sas[0][:nAvoid],c='gray',alpha=0.3,label='Samples')
    for i,sample in tqdm(enumerate(folders), desc='Plotting simulation data'):
        axes[0].plot(x[:nAvoid],Scs[i][:nAvoid],c='gray',alpha=0.3)
        axes[1].plot(x[:nAvoid],Sbs[i][:nAvoid],c='gray',alpha=0.3)
        axes[2].plot(x[:nAvoid],Sas[i][:nAvoid],c='gray',alpha=0.3)
    
    # CT scan profile
    if isCTscan:
    #    axes[0].scatter(CTscan_results['x_SoCT'],CTscan_results['So_CT'],c='k',label='CT scan') 
    #    axes[1].scatter(CTscan_results['x_SwCT'],CTscan_results['Sw_CT'],c='k',label='CT scan')
    #    axes[2].scatter(CTscan_results['x_SgCT'],CTscan_results['Sg_CT'],c='k',label='CT scan')
        axes[0].plot(CTscan_results['x_SoCT'],CTscan_results['So_CT'],'k.-', markersize=15,label='CT scan') 
        axes[1].plot(CTscan_results['x_SwCT'],CTscan_results['Sw_CT'],'k.-', markersize=15,label='CT scan')
        axes[2].plot(CTscan_results['x_SgCT'],CTscan_results['Sg_CT'],'k.-', markersize=15,label='CT scan')
    
    if isDARTS:
        axes[0].plot(DARTS_results['x_So_DARTS'],DARTS_results['So_DARTS'],'r.-',label='DARTS')
        axes[1].plot(DARTS_results['x_Sw_DARTS'],DARTS_results['Sw_DARTS'],'b.-',label='DARTS')
        axes[2].plot(DARTS_results['x_Sg_DARTS'],DARTS_results['Sg_DARTS'],'g.-',label='DARTS')

    axes[0].set_xlabel('x [m]')
    axes[0].set_ylabel(r'$S_o$ [-]')
    axes[0].set_ylim([0,1])
    axes[0].grid()
    axes[0].legend()
    axes[1].set_xlabel('x [m]')
    axes[1].set_ylabel(r'$S_w$ [-]')
    axes[1].set_ylim([0,1])
    axes[1].grid()
    axes[1].legend()
    axes[2].set_xlabel('x [m]')
    axes[2].set_ylabel(r'$S_g$ [-]')
    axes[2].set_ylim([0,1])
    axes[2].grid()
    axes[2].legend()
    plt.tight_layout()
    plt.savefig(experiment_dir + '/profiles_PVI_'+str(tD)+'.png',dpi=200)
    plt.show()
    plt.close()


    # PLOT F STARS:
    Fdrys = []
    for i in range(len(X_ED)):
        Fdrys.append(Fdry(Sbs[i][:], X_ED.iloc[i,:].to_numpy()))

    SF_values = X_ED.iloc[:,4].to_numpy()
    norm = Normalize(vmin=np.min(SF_values), vmax=np.max(SF_values))
    cmap = cm.get_cmap('viridis')

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i,sample in tqdm(enumerate(folders), desc='Plotting simulation data'):
        color = cmap(norm(SF_values[i]))
        axes[0].plot(x[1:nAvoid],Foils[i][1:nAvoid],color=color,alpha=0.3)
        axes[1].plot(x[1:nAvoid],Fshears[i][1:nAvoid],color=color,alpha=0.3)
        axes[2].plot(x[1:nAvoid],Fdrys[i][1:nAvoid],color=color,alpha=0.3)

    axes[0].set_xlabel('x [m]')
    axes[0].set_ylabel(r'$F_{oil}$ [-]')
    axes[0].set_ylim([-0.05,1.05])
    axes[0].grid()
    axes[0].legend()
    axes[1].set_xlabel('x [m]')
    axes[1].set_ylabel(r'$F_{shear}$ [-]')
    axes[1].set_ylim([-0.05,1.05])
    axes[1].grid()
    axes[1].legend()
    axes[2].set_xlabel('x [m]')
    axes[2].set_ylabel(r'$F_{dry}$ [-]')
    axes[2].set_ylim([-0.05,1.05])
    axes[2].grid()
    axes[2].legend()

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=axes[2])
    cbar.set_label(r"$S_w^{\ast}$")

    plt.tight_layout()
    plt.savefig(experiment_dir + '/FStars_profiles_PVI_'+str(tD)+'.png',dpi=200)
    plt.show()
    plt.close()

    # PLOT SATURATIONS
    SF_values = X_ED.iloc[:,4].to_numpy()
    norm = Normalize(vmin=np.min(SF_values), vmax=np.max(SF_values))
    cmap = cm.get_cmap('viridis')
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    for i,sample in tqdm(enumerate(folders), desc='Plotting simulation data'):
        color = cmap(norm(SF_values[i]))
        ax.plot(x[:nAvoid],Scs[i][:nAvoid],color=color)
    if isDARTS:
        ax.plot(DARTS_results['x_So_DARTS'],DARTS_results['So_DARTS'],'r.-',label='DARTS')
    if isCTscan:
        ax.plot(CTscan_results['x_SoCT'],CTscan_results['So_CT'],'k.-', markersize=15,label='CT scan') 
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, ax=ax)
    # cbar.set_label(r"$fmoil$")
    plt.xlabel('x [m]')
    plt.ylabel(r'$S_{o}$ [-]')
    plt.ylim([-0.05,1.05])
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(experiment_dir + '/So_fmoil_profiles_PVI_'+str(tD)+'.png',dpi=200)
    plt.show()


    # # PLOT HISTOGRAM OF MEDIA SATURATIONS
    # Sw_median = []
    # for i,sample in tqdm(enumerate(folders), desc='Plotting simulation data'):
    #     Sw_median.append(np.median(Sbs[i][:nAvoid]))
    
    # print(Sw_median)

    # plt.figure(figsize=(5,5))
    # plt.hist(np.array(Sw_median),bins=10)
    # plt.xlabel(r"$S_w$")
    # plt.ylabel(r"Frequency")
    # plt.grid()
    # plt.tight_layout()
    # plt.show()