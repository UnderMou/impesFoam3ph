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

    # Inputs:
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_GammaBeta'

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
    tf = 40000
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
    # folders = folders[:100]
    # indexes_to_remove = [49,476,733,953,914,1443,1451]    # original
    # indexes_to_remove = [18, 33, 49, 57, 106, 129, 177, 182, 198, 222, 227, 236, 245, 246, 284, 305, 310, 388, 402, 412, 416, 444, 456, 460, 476, 482, 493, 494, 518, 525, 580, 603, 625, 649, 670, 678, 698, 721, 726, 732, 783, 803, 835, 849, 896, 913, 915, 942, 952, 962, 998, 1015, 1037, 1086, 1095, 1098, 1102, 1125, 1136, 1147, 1171, 1193, 1219, 1275, 1286]
    indexes_to_remove = [0, 17, 18, 49, 57, 71, 78, 84, 94, 106, 125, 128, 129, 137, 142, 159, 162, 166, 175, 177, 182, 198, 201, 222, 227, 236, 245, 246, 255, 259, 279, 284, 286, 289, 305, 310, 312, 345, 364, 377, 383, 388, 401, 402, 404, 406, 412, 416, 426, 444, 456, 460, 463, 468, 476, 482, 493, 494, 498, 504, 518, 525, 527, 545, 560, 571, 573, 580, 603, 614, 615, 625, 643, 647, 649, 662, 670, 678, 690, 695, 696, 697, 698, 705, 706, 721, 726, 732, 738, 770, 772, 783, 792, 795, 803, 817, 822, 835, 849, 867, 869, 896, 904, 911, 912, 913, 915, 922, 942, 952, 962, 970, 972, 983, 998, 1015, 1037, 1067, 1070, 1086, 1095, 1098, 1102, 1119, 1125, 1136, 1142, 1147, 1155, 1171, 1172, 1186, 1193, 1194, 1203, 1209, 1214, 1216, 1219, 1224, 1228, 1238, 1244, 1245, 1275, 1286]
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


    # Check oscilations
    nAvoid = -1
    indexes_to_remove = []
    for i,sample in enumerate(folders):
        count = checkOscilations(prod_w[i][:nAvoid])
        if count >= 5:
            indexes_to_remove.append(i)
        print(i, count)

    print(indexes_to_remove)
