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
        t_find = 750 
        isCTscan = True
        isDARTS = True
    elif tD in [.82, 1.0]:  
        t_find = 1705
        isCTscan = True
        isDARTS = True
    elif tD == 10.40: 
        raise KeyError("The 'tD' provided is not expected")
        t_find = 4975
        isCTscan = True
        isDARTS = False
    else: raise KeyError("The 'tD' provided is not expected")
    return t_find, isCTscan, isDARTS

if __name__ == '__main__':

    # tDs = [0.36, 0.44]
    tDs = [0.82, 1.0]
    isDARTS = [True, False]
    isCTscan = [True, True]
    colors = ['g', 'r']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i,tD in enumerate(tDs):
        if isCTscan[i]:
            Sg_CT_file = 'Tang/Sg_PVI_'+str(tD)+'.csv'
            Sw_CT_file = 'Tang/Sw_PVI_'+str(tD)+'.csv'
            So_CT_file = 'Tang/So_PVI_'+str(tD)+'.csv'
            CTscan_results = get_CTscanData(So_CT_file, Sg_CT_file, Sw_CT_file)

            axes[0].plot(CTscan_results['x_SoCT'],CTscan_results['So_CT'],f'{colors[i]}.-', markersize=15,label=f'CT scan {tD}') 
            axes[1].plot(CTscan_results['x_SwCT'],CTscan_results['Sw_CT'],f'{colors[i]}.-', markersize=15,label=f'CT scan {tD}')
            axes[2].plot(CTscan_results['x_SgCT'],CTscan_results['Sg_CT'],f'{colors[i]}.-', markersize=15,label=f'CT scan {tD}')

        if isDARTS[i]:
            Sg_DARTS_file = 'Tang/Sg_PVI_'+str(tD)+'_DARTS.csv'
            Sw_DARTS_file = 'Tang/Sw_PVI_'+str(tD)+'_DARTS.csv'
            So_DARTS_file = 'Tang/So_PVI_'+str(tD)+'_DARTS.csv'
            DARTS_results = get_DARTSdata(Sg_DARTS_file, Sw_DARTS_file, So_DARTS_file)

            axes[0].plot(DARTS_results['x_So_DARTS'],DARTS_results['So_DARTS'],f'{colors[i]}.--',label=f'DARTS {tD}')
            axes[1].plot(DARTS_results['x_Sw_DARTS'],DARTS_results['Sw_DARTS'],f'{colors[i]}.--',label=f'DARTS {tD}')
            axes[2].plot(DARTS_results['x_Sg_DARTS'],DARTS_results['Sg_DARTS'],f'{colors[i]}.--',label=f'DARTS {tD}')

    axes[0].set_xlabel('x [m]')
    axes[0].set_ylabel(r'$S_o$ [-]')
    axes[0].grid()
    axes[0].legend()
    axes[1].set_xlabel('x [m]')
    axes[1].set_ylabel(r'$S_w$ [-]')
    axes[1].grid()
    axes[1].legend()
    axes[2].set_xlabel('x [m]')
    axes[2].set_ylabel(r'$S_g$ [-]')
    axes[2].grid()
    axes[2].legend()
    plt.tight_layout()
    plt.show()