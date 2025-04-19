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

from aux_funcs import *

def checkOscilations(S):
    L2 = int(len(S)/2)
    
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

    # Inputs:
    # experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_fitted_v2'
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_fitted_fixedSF'

    # Reservoir:
    reservoir = {
        'H': 2 * 0.01772453850905516,
        'L': 0.4,
        'phi': 0.22
    }
    reservoir['A'] = reservoir['H']**2
    reservoir['PoreVol'] = reservoir['phi'] * (reservoir['A']*reservoir['L'])

    # Time controls:
    ti = 0
    tf = 10000
    write_interval = 200
    t = np.linspace(ti,tf,int(tf/write_interval + 1))

    # Residual saturations:
    residuals = {
        'Swc': 0,#0.197,
        'Sgr': 0,#0.013,
        'Sor': 0 #0.103
    }

    # Sampling results
    experiment_dir = set_experimentDir(experiment_name)
    folders = get_folders(experiment_dir)
    folders = folders[:500]
    # indexes_to_remove = [55, 63, 65, 66, 136, 137, 175, 184, 205, 226, 232, 247, 311, 381, 384, 433]    # original
    indexes_to_remove = [18, 19, 21, 36, 49, 57, 66, 182, 222, 236, 245, 284, 305, 310, 340, 388, 402, 412, 416, 456, 482, 493, 494] # fixedSF
    folders = remove_indexes(folders, indexes_to_remove) 
    nsamples = len(folders)

    
    # EVALUATING QoIs:
    nAvoid = -2
    time_steps = len(t[1:])

    # 1st QoI - Cumulative phase production  
    prod_w, prod_o, prod_g = eval_production(experiment_dir, 
                                             folders,
                                             reservoir,
                                             nsamples,
                                             time_steps,
                                             nAvoid,
                                             residuals,
                                             write_interval)
    plot_production(nsamples, prod_o, prod_w, prod_g, t, experiment_dir)

    cuml_prod_w, cuml_prod_o, cuml_prod_g = eval_cumulativeProduction(prod_w,
                                                                      prod_o,
                                                                      prod_g,
                                                                      nsamples,
                                                                      time_steps)
    plot_Cumulativeproduction(nsamples, cuml_prod_o, cuml_prod_w, cuml_prod_g, t, experiment_dir)                                                                      



    # 2nd QoI - Pressure drop
    pressDrop = eval_pressureDrop(experiment_dir,
                                  folders,
                                  nsamples,
                                  time_steps,
                                  nAvoid)
    plot_pressDrop(nsamples, pressDrop, t, experiment_dir)



    # 3rd QoI - Oil recovery factor 
    Soi = 0.46
    OOIP = Soi
    # print(OOIP)
    ORF = (cuml_prod_o / OOIP) * 100
    plot_OilRecFactor(nsamples, ORF, t, experiment_dir)


    # 4th QoI - Oil cut
    totalProd = prod_w + prod_o + prod_g
    Oil_cut = (prod_o / totalProd) * 100
    plot_OilCut(nsamples, Oil_cut, t, experiment_dir)

    # Plot QoI 2, 3 and 4
    plot_PD_ORF_OC(nsamples, t, experiment_dir, pressDrop, ORF, Oil_cut)

    # Plot histograms
    idt = -1 # steady state index
    plot_QoIs_histograms_steadyState(idt,
                                     experiment_dir, 
                                     cuml_prod_o,
                                     cuml_prod_w,
                                     cuml_prod_g,
                                     pressDrop,
                                     ORF,
                                     nBins = 30)
    
    # Saving QoIs
    # 1st QoI - Cumulative phase production
    np.savetxt(experiment_dir + '/QoI_1_prodOil.csv', prod_o, delimiter=',')
    np.savetxt(experiment_dir + '/QoI_1_prodWater.csv', prod_w, delimiter=',')
    np.savetxt(experiment_dir + '/QoI_1_prodGas.csv', prod_g, delimiter=',')
    np.savetxt(experiment_dir + '/QoI_1_cumulativeProdOil.csv', cuml_prod_o, delimiter=',')
    np.savetxt(experiment_dir + '/QoI_1_cumulativeProdWater.csv', cuml_prod_w, delimiter=',')
    np.savetxt(experiment_dir + '/QoI_1_cumulativeProdGas.csv', cuml_prod_g, delimiter=',')

    # 2nd QoI - Pressure drop
    np.savetxt(experiment_dir + '/QoI_2_pressDrop.csv', pressDrop, delimiter=',')

    # 3rd QoI - Oil recovery factor
    np.savetxt(experiment_dir + '/QoI_3_OilRecFac.csv', ORF, delimiter=',')

    # 4th QoI - Oil cut
    np.savetxt(experiment_dir + '/QoI_4_OilCut.csv', Oil_cut, delimiter=',')


    # Check oscilations
    # nAvoid = -1
    indexes_to_remove = []
    for i,sample in enumerate(folders):
        # count = checkOscilations(pressDrop[i][:nAvoid])
        count = checkOscilations(pressDrop[i][:nAvoid])
        if count > 2:
            indexes_to_remove.append(i)
        print(i, count)

    print(indexes_to_remove)
