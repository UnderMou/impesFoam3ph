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

from uqpylab import sessions, display_util

from aux_funcs import *

def tFind(tD):
    if tD in [.36]:  
        t_find = 800 
        isCTscan = True
        isDARTS = True
    elif tD in [.82]:  
        t_find = 1700
        isCTscan = True
        isDARTS = True
    elif tD in [.6]:
        t_find = 1400 
        isCTscan = False
        isDARTS = False
    elif tD in [.6]:
        t_find = 1400 
        isCTscan = False
        isDARTS = False
    else: 
        raise KeyError("The 'tD' provided is not expected")
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

    # Inputs
    experiment_name = 'SA_CTscan_ViscoFoil'
    nSamples = 1000
    ti = 0
    tf = 2000
    write_interval = 200
    t = np.linspace(ti,tf,int(tf/write_interval + 1))

    # Get samples from OpenFOAM simulations parameters
    experiment_dir = set_experimentDir(experiment_name)
    X_ED = pd.read_csv(experiment_dir + '/X_ED.csv',header=None)
    indexes_to_remove = [32, 46, 64, 109, 166, 197, 258, 318, 376, 736, 858]
    X_ED = X_ED.drop(indexes_to_remove).reset_index(drop=True)
    X_ED = X_ED[:nSamples]
    print(X_ED.shape)
    print(X_ED.head())

    # # #########################################
    # # EVAL OIL BANK QoI FROM OIL CUT/PRODUCTION
    # # #########################################

    # Prod_o = pd.read_csv(experiment_dir + '/QoI_1_prodOil.csv', header=None).iloc[:nSamples,:]
    # Prod_w = pd.read_csv(experiment_dir + '/QoI_1_prodWater.csv', header=None).iloc[:nSamples,:]
    # Prod_g = pd.read_csv(experiment_dir + '/QoI_1_prodGas.csv', header=None).iloc[:nSamples,:]
    # # QoI_2_SIM_PressDrop = pd.read_csv(experiment_dir + '/QoI_2_pressDrop.csv', header=None).iloc[:nSamples,:]
    # # QoI_3_SIM_OilRecFac = pd.read_csv(experiment_dir + '/QoI_3_OilRecFac.csv', header=None).iloc[:nSamples,:]
    # # OilCut = pd.read_csv(experiment_dir + '/QoI_4_OilCut.csv', header=None).iloc[:nSamples,:]
    
    # totalProd = Prod_w + Prod_o + Prod_g
    # # OilCut = (Prod_o / totalProd) * 100 # [%]
    # OilCut = Prod_o # [PV]

    # OilBank_Area = eval_oilBank_AreaBelow_fromOilCut(OilCut, t, write_interval)
    # OilBank_PulseDur = eval_oilBank_PulseDuration_fromOilCut(OilCut, t)
    # OilBank_PulseHei = eval_oilBank_PulseHeight_fromOilCut(OilCut, t)
    # plot_oilBank(OilBank_PulseDur, OilBank_PulseHei, OilBank_Area)


    # ##############################################
    # EVAL OIL BANK QoI FROM OIL SATURATION PROFILES
    # ##############################################

    tD = 0.36
    # indexes_to_remove = [49,476,733,953,914,1443,1451]    # original
    nAvoid = -1

    t_find, isCTscan, isDARTS = tFind(tD)
    idx = np.argmin(np.abs(t-t_find))

    # Residual saturations:
    residuals = {
        'Swc': 0.197,
        'Sgr': 0.013,
        'Sor': 0.103
    }

    Scs = np.loadtxt(experiment_dir + '/Scs_' + str(tD) + '.csv')
    nCells = len(Scs[0])
    x = np.linspace(0,1,nCells)
    L = x[-1]
    x = x/L

    OilBank_length, init_idxs, end_idxs = eval_oilBank_length_fromOilSat(pd.DataFrame(Scs), x, Soi = 0.46, nAvoid=nAvoid)
    OilBank_height = eval_oilBank_height_fromOilSat(pd.DataFrame(Scs), x, init_idxs, end_idxs, Soi = 0.46)
    OilBank_area = eval_oilBank_AreaBelow_fromOilCut(pd.DataFrame(Scs), x, init_idxs, end_idxs, Soi = 0.46)


    np.savetxt(experiment_dir + '/OilBank_length_OilSat_' + str(tD) + '.csv', OilBank_length)
    np.savetxt(experiment_dir + '/OilBank_height_OilSat_' + str(tD) + '.csv', OilBank_height)
    np.savetxt(experiment_dir + '/OilBank_area_OilSat_' + str(tD) + '.csv', OilBank_area)

    fig, axs = plt.subplots(1,3, figsize=(10,5))

    axs[0].hist(OilBank_length, bins=7)
    axs[0].set_xlabel('Oil bank length [xD = x/L]')
    axs[0].set_ylabel('Frequency')
    axs[0].grid()

    axs[1].hist(OilBank_height, bins=7)
    axs[1].set_xlabel(r'Oil bank Height [$S_o$]')
    axs[1].set_ylabel('Frequency')
    axs[1].grid()

    axs[2].hist(OilBank_area, bins=7)
    axs[2].set_xlabel(r'Oil bank area [$S_o$]')
    axs[2].set_ylabel('Frequency')
    axs[2].grid()

    plt.tight_layout()
    plt.savefig(experiment_dir + '/QoIs_histogram_oilBank_' + str(tD) + '.png', dpi=200)
    plt.show()

    plot_fromOilSat(pd.DataFrame(Scs), x, init_idxs, end_idxs, Soi = 0.46)
