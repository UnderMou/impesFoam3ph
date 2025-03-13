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

from uqpylab import sessions, display_util

from aux_funcs import *

if __name__ == '__main__':

    # Inputs
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples'
    nSamples = 736
    ti = 0
    tf = 40000
    write_interval = 200
    nPCEs = 1
    trainingSet_size = 0.75
    times2PCE = [39800] # np.linspace(2000,tf,nPCEs).astype(int)

    # Get samples from OpenFOAM simulations parameters
    experiment_dir = set_experimentDir(experiment_name)
    X_ED = pd.read_csv(experiment_dir + '/X_ED.csv',header=None)
    X_ED = X_ED[:nSamples]
    
    index = round(X_ED.shape[0]*trainingSet_size)
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    idxs_times2PCE = [np.argmin(np.abs(t-time)) for time in times2PCE]
    print(idxs_times2PCE)

    # Initialize QoIs array from simulations
    QoI_1_SIM_CumProd_o = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdOil.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_1_SIM_CumProd_w = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdWater.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_1_SIM_CumProd_g = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdGas.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_2_SIM_PressDrop = pd.read_csv(experiment_dir + '/QoI_2_pressDrop.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_3_SIM_OilRecFac = pd.read_csv(experiment_dir + '/QoI_3_OilRecFac.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_4_SIM_OilCut = pd.read_csv(experiment_dir + '/QoI_4_OilCut.csv', header=None).iloc[:nSamples,idxs_times2PCE]

    # Initialize PCEs list
    QoI_1_PCE_CumProd_o = []
    QoI_1_PCE_CumProd_w = []
    QoI_1_PCE_CumProd_g = []
    QoI_2_PCE_PressDrop = []
    QoI_3_PCE_OilRecFac = []
    QoI_4_PCE_OilCut = []

    # Initialize sensitivity analysis Sobol indexes list
    QoI_1_SA_CumProd_o = []
    QoI_1_SA_CumProd_w = []
    QoI_1_SA_CumProd_g = []
    QoI_2_SA_PressDrop = []
    QoI_3_SA_OilRecFac = []
    QoI_4_SA_OilCut = []

    # Initialize common plotting parameters
    display_util.load_plt_defaults()
    uq_colors = display_util.get_uq_color_order()

    # Start the session
    myToken = '270fe799b4e4e4f004da79cec3dffb07cf9da78a' # The user's token to access the UQCloud API
    UQCloud_instance = 'https://uqcloud.ethz.ch' # The UQCloud instance to use
    mySession = sessions.cloud(host=UQCloud_instance, token=myToken)
    # (Optional) Get a convenient handle to the command line interface
    uq = mySession.cli
    # Reset the session
    mySession.reset()

    # Set random seed for reproducibility
    uq.rng(0,'twister')

    # Reference values
    fmmobR = 50000.0
    sfdryR = 0.215
    sfbetR = 19950.0
    fmcapR = 1e-4
    epcapR = 1.321
    fmoilR = 0.823
    delta = 0.8
    epoilR = 3.827

    # Set parameters bounds
    dpr = 0.2
    dists = {
    'num_vars': 8,
    'names': ['fmmob', 'sfdry', 'sfbet', 'fmcap', 'epcap', 'fmoil', 'delta', 'epoil'],
    'bounds': [[fmmobR*(1-dpr), fmmobR*(1+dpr)],
               [sfdryR*(1-dpr), sfdryR*(1+dpr)],
               [sfbetR*(1-dpr), sfbetR*(1+dpr)],
               [fmcapR*(1-dpr), fmcapR*(1+dpr)],
               [epcapR*(1-dpr), epcapR*(1+dpr)],
               [fmoilR*(1-dpr), fmoilR*(1+dpr)],
               [delta*(1-dpr), delta*(1+dpr)],
               [epoilR*(1-dpr), epoilR*(1+dpr)],
              ]
    }

    # Set marginals
    InputOpts = {
        "Marginals": [
            {"Name": dists['names'][0],
            "Type": "Uniform",
            "Parameters": dists['bounds'][0]
            },
            {"Name": dists['names'][1],
            "Type": "Uniform",
            "Parameters": dists['bounds'][1]
            },
            {"Name": dists['names'][2],
            "Type": "Uniform",
            "Parameters": dists['bounds'][2]
            },
            {"Name": dists['names'][3],
            "Type": "Uniform",
            "Parameters": dists['bounds'][3]
            },
            {"Name": dists['names'][4],
            "Type": "Uniform",
            "Parameters": dists['bounds'][4]
            },
            {"Name": dists['names'][5],
            "Type": "Uniform",
            "Parameters": dists['bounds'][5]
            },
            {"Name": dists['names'][6],
            "Type": "Uniform",
            "Parameters": dists['bounds'][6]
            },
            {"Name": dists['names'][7],
            "Type": "Uniform",
            "Parameters": dists['bounds'][7]
            }
        ]
    }

    myInput = uq.createInput(InputOpts)
    uq.print(myInput)

    # Set sensitivity analysis
    PCESobol = {
        'Type': 'Sensitivity',
        'Method': 'Sobol',
        'Sobol': {
            'Order': 5
        }
    }

    # Training PCEs
    for i in range(nPCEs):
        print("PCE model no." + str(i+1) + "/" + str(nPCEs))
        print("Time: " + str(t[idxs_times2PCE[i]]))

        # QoI 1 - Cumulative oil production
        print('QoI 1 - Cumulative oil production')
        QoI_1_PCE_CumProd_o.append(
            eval_PCE(uq, 
                     X = X_ED, 
                     Y = QoI_1_SIM_CumProd_o,
                     training_idx = index, 
                     QoI_col = i)
        )
        uq.print(QoI_1_PCE_CumProd_o[-1])
        check_pceAcc(uq,
                     pce_model = QoI_1_PCE_CumProd_o[-1],
                     X_SIM = X_ED.iloc[index:, :].to_numpy().tolist(),
                     Y_SIM = QoI_1_SIM_CumProd_o.iloc[index:,i],
                     saveDir=experiment_dir + '/PCEsGraphs',
                     fileName = 'QoI_1_CumProd_o.png'
        )
        QoI_1_SA_CumProd_o.append(uq.createAnalysis(PCESobol))
        print(QoI_1_SA_CumProd_o)
        print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_1_SA_CumProd_o, PCESobol, idx = -1)}')

        # QoI 1 - Cumulative water production
        print('QoI 1 - Cumulative water production')
        QoI_1_PCE_CumProd_w.append(
            eval_PCE(uq, 
                     X = X_ED, 
                     Y = QoI_1_SIM_CumProd_w,
                     training_idx = index, 
                     QoI_col = i)
        )
        uq.print(QoI_1_PCE_CumProd_w[-1])
        check_pceAcc(uq,
                     pce_model = QoI_1_PCE_CumProd_w[-1],
                     X_SIM = X_ED.iloc[index:, :].to_numpy().tolist(),
                     Y_SIM = QoI_1_SIM_CumProd_w.iloc[index:,i],
                     saveDir=experiment_dir + '/PCEsGraphs',
                     fileName = 'QoI_1_CumProd_w.png'
        )
        QoI_1_SA_CumProd_w.append(uq.createAnalysis(PCESobol))
        print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_1_SA_CumProd_w, PCESobol, idx = -1)}')

        # QoI 1 - Cumulative gas production
        print('QoI 1 - Cumulative gas production')
        QoI_1_PCE_CumProd_g.append(
            eval_PCE(uq, 
                     X = X_ED, 
                     Y = QoI_1_SIM_CumProd_g,
                     training_idx = index, 
                     QoI_col = i)
        )
        uq.print(QoI_1_PCE_CumProd_g[-1])
        check_pceAcc(uq,
                     pce_model = QoI_1_PCE_CumProd_g[-1],
                     X_SIM = X_ED.iloc[index:, :].to_numpy().tolist(),
                     Y_SIM = QoI_1_SIM_CumProd_g.iloc[index:,i],
                     saveDir=experiment_dir + '/PCEsGraphs',
                     fileName = 'QoI_1_CumProd_g.png'
        )
        QoI_1_SA_CumProd_g.append(uq.createAnalysis(PCESobol))
        print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_1_SA_CumProd_g, PCESobol, idx = -1)}')

        # QoI 2 - Pressure drop
        print('QoI 2 - Pressure drop')
        QoI_2_PCE_PressDrop.append(
            eval_PCE(uq, 
                     X = X_ED, 
                     Y = QoI_2_SIM_PressDrop,
                     training_idx = index, 
                     QoI_col = i)
        )
        uq.print(QoI_2_PCE_PressDrop[-1])
        check_pceAcc(uq,
                     pce_model = QoI_2_PCE_PressDrop[-1],
                     X_SIM = X_ED.iloc[index:, :].to_numpy().tolist(),
                     Y_SIM = QoI_2_SIM_PressDrop.iloc[index:,i],
                     saveDir=experiment_dir + '/PCEsGraphs',
                     fileName = 'QoI_2_PressDrop.png'
        )
        QoI_2_SA_PressDrop.append(uq.createAnalysis(PCESobol))
        print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_2_SA_PressDrop, PCESobol, idx = -1)}')

        # QoI 3 - Oil recovery factor
        print('QoI 3 - Oil recovery factor')
        QoI_3_PCE_OilRecFac.append(
            eval_PCE(uq, 
                     X = X_ED, 
                     Y = QoI_3_SIM_OilRecFac,
                     training_idx = index, 
                     QoI_col = i)
        )
        uq.print(QoI_3_PCE_OilRecFac[-1])
        check_pceAcc(uq,
                     pce_model = QoI_3_PCE_OilRecFac[-1],
                     X_SIM = X_ED.iloc[index:, :].to_numpy().tolist(),
                     Y_SIM = QoI_3_SIM_OilRecFac.iloc[index:,i],
                     saveDir=experiment_dir + '/PCEsGraphs',
                     fileName = 'QoI_3_OilRecFac.png'
        )
        QoI_3_SA_OilRecFac.append(uq.createAnalysis(PCESobol))
        print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_3_SA_OilRecFac, PCESobol, idx = -1)}')

    # Plot Sobol indixes
    nOrder = PCESobol['Sobol']['Order']
    plot_SobolIndex(nOrder, QoI_1_SA_CumProd_o, saveDir=experiment_dir + '/SA_analysis/QoI_1_CumProd_o')
    plot_SobolIndex(nOrder, QoI_1_SA_CumProd_w, saveDir=experiment_dir + '/SA_analysis/QoI_1_CumProd_w')
    plot_SobolIndex(nOrder, QoI_1_SA_CumProd_g, saveDir=experiment_dir + '/SA_analysis/QoI_1_CumProd_g')
    plot_SobolIndex(nOrder, QoI_2_SA_PressDrop, saveDir=experiment_dir + '/SA_analysis/QoI_2_SA_PressDrop')
    plot_SobolIndex(nOrder, QoI_3_SA_OilRecFac, saveDir=experiment_dir + '/SA_analysis/QoI_3_SA_OilRecFac')
