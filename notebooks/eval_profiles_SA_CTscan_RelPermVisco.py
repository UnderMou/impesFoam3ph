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

from aux_funcs import *

def tFind(tD):
    if tD in [.36, .44]:  
        t_find = 800 #750
        isCTscan = True
        isDARTS = True
    elif tD in [.82, 1.0]:  
        t_find = 1700
        isCTscan = True
        isDARTS = True
    elif tD in [.6]:
        t_find = 1400 
        isCTscan = False
        isDARTS = False
    elif tD in [.2]:
        t_find = 400 
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

    tD = 0.36
    experiment_name = 'SA_CTscan_RelPermVisco'

    # Residual saturations:
    residuals = {
        'Swc': 0.197,
        'Sgr': 0.013,
        'Sor': 0.103
    }

    # Time controls:
    ti = 0
    tf = 2000 
    write_interval = 200
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    t_find, isCTscan, isDARTS = tFind(tD)
    # t_find, isCTscan, isDARTS = 600, False, False
    idx = np.argmin(np.abs(t-t_find))
    print(t_find)

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
    # folders = folders[:251]
    indexes_to_remove = []
    folders = remove_indexes(folders, indexes_to_remove) 
    Sas, Sbs, Scs = get_Saturations(folders, experiment_dir, idx, residuals)
    X_ED = get_samples(experiment_dir)
    nCells = len(Sas[0])
    x = np.linspace(0,1,nCells)
    L = x[-1]
    x = x/L

    # # Check oscilations
    # nAvoid = -1
    # indexes_to_remove = []
    # for i,sample in enumerate(folders):
    #     count = checkOscilations(Scs[i][:nAvoid])
    #     if count >= 5:
    #         indexes_to_remove.append(i)
    #     print(i, count)

    print(indexes_to_remove)
    folders = remove_indexes(folders, indexes_to_remove) 
    Sas = remove_indexes(Sas, indexes_to_remove) 
    Sbs = remove_indexes(Sbs, indexes_to_remove) 
    Scs = remove_indexes(Scs, indexes_to_remove) 

    np.savetxt(experiment_dir + '/Scs_' + str(tD) + '.csv', Scs)

    # Plotting:

    # filtering
    idX = np.argmin(np.abs(x-0.373))
    idY = np.argmin(np.abs(x-0.123))
    X_ED = pd.DataFrame(X_ED)

    nAvoid = -1
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    # Samples profile
    for i,sample in tqdm(enumerate(folders), desc='Plotting simulation data'):
        # axes[0].plot(x[:nAvoid],Scs[i][:nAvoid])
        # axes[1].plot(x[:nAvoid],Sbs[i][:nAvoid])
        # axes[2].plot(x[:nAvoid],Sas[i][:nAvoid])

        if Scs[i][idX] > 0.8 and Scs[i][idY] > 0.21:
            print(i, X_ED.iloc[i,:])
            axes[0].plot(x[:nAvoid],Scs[i][:nAvoid])
            axes[1].plot(x[:nAvoid],Sbs[i][:nAvoid])
            axes[2].plot(x[:nAvoid],Sas[i][:nAvoid])

        # axes[0].plot(x[:nAvoid],Scs[i][:nAvoid],c='gray',alpha=0.3)
        # axes[1].plot(x[:nAvoid],Sbs[i][:nAvoid],c='gray',alpha=0.3)
        # axes[2].plot(x[:nAvoid],Sas[i][:nAvoid],c='gray',alpha=0.3)
    
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