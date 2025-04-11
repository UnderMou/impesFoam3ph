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

if __name__ == '__main__':

    # Inputs:
    experiment_name = 'SA_CTscan_RelPermVisco'

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
    tf = 2000
    write_interval = 200
    t = np.linspace(ti,tf,int(tf/write_interval + 1))

    # Residual saturations:
    residuals = {
        'Swc': 0.197,
        'Sgr': 0.013,
        'Sor': 0.103
    }

    # Sampling results
    experiment_dir = set_experimentDir(experiment_name)
    folders = get_folders(experiment_dir)
    folders = folders[:257]
    # indexes_to_remove = [49,476,733,953,914,1443,1451]    # original
    indexes_to_remove = []
    folders = remove_indexes(folders, indexes_to_remove) 
    nsamples = len(folders)
    
    # EVALUATING QoIs:
    nAvoid = -1
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
    Soi = 0.41902
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
