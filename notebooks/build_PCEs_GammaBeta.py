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

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams.update({'font.size': 16})

    # Inputs
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_GammaBeta'
    nSamples = 1303
    indexes_to_remove = [0, 17, 18, 49, 57, 71, 78, 84, 94, 106, 125, 128, 129, 137, 142, 159, 162, 166, 175, 177, 182, 198, 201, 222, 227, 236, 245, 246, 255, 259, 279, 284, 286, 289, 305, 310, 312, 345, 364, 377, 383, 388, 401, 402, 404, 406, 412, 416, 426, 444, 456, 460, 463, 468, 476, 482, 493, 494, 498, 504, 518, 525, 527, 545, 560, 571, 573, 580, 603, 614, 615, 625, 643, 647, 649, 662, 670, 678, 690, 695, 696, 697, 698, 705, 706, 721, 726, 732, 738, 770, 772, 783, 792, 795, 803, 817, 822, 835, 849, 867, 869, 896, 904, 911, 912, 913, 915, 922, 942, 952, 962, 970, 972, 983, 998, 1015, 1037, 1067, 1070, 1086, 1095, 1098, 1102, 1119, 1125, 1136, 1142, 1147, 1155, 1171, 1172, 1186, 1193, 1194, 1203, 1209, 1214, 1216, 1219, 1224, 1228, 1238, 1244, 1245, 1275, 1286]
    nSamples -= len(indexes_to_remove)
    ti = 0
    tf = 40000
    write_interval = 200
    nPCEs = 1
    trainingSet_size = 0.75
    times2PCE = [39800] # np.linspace(2000,tf,nPCEs).astype(int)

    # Get samples from OpenFOAM simulations parameters
    experiment_dir = set_experimentDir(experiment_name)
    X_ED = pd.read_csv(experiment_dir + '/X_ED.csv',header=None)
    # indexes_to_remove = [49,476,733,953,914,1443,1451]    original
    # indexes_to_remove = [18, 33, 49, 57, 106, 129, 177, 182, 198, 222, 227, 236, 245, 246, 284, 305, 310, 388, 402, 412, 416, 444, 456, 460, 476, 482, 493, 494, 518, 525, 580, 603, 625, 649, 670, 678, 698, 721, 726, 732, 783, 803, 835, 849, 896, 913, 915, 942, 952, 962, 998, 1015, 1037, 1086, 1095, 1098, 1102, 1125, 1136, 1147, 1171, 1193, 1219, 1275, 1286]
    X_ED = X_ED.drop(indexes_to_remove).reset_index(drop=True)
    X_ED = X_ED[:nSamples]
    print(X_ED.shape)
    print(X_ED.head())
    # X_ED.iloc[:, 5] = -(X_ED.iloc[:,5]/X_ED.iloc[:,4] - 1)
    # floil_samples = fmoil_samples * (1.0 - delta_samples)
    print(X_ED.head())

    index = round(X_ED.shape[0]*trainingSet_size)
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    idxs_times2PCE = [np.argmin(np.abs(t-time)) for time in times2PCE]
    print(idxs_times2PCE)

    # Initialize QoIs array from simulations
    QoI_1_SIM_CumProd_o = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdOil.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    print(QoI_1_SIM_CumProd_o.shape)
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
    names = ['fmmob', 'sfdry', 'sfbet', 'epcap', 'fmoil', 'floil', 'epoil']
    mu_fmmob  = 50000.0
    mu_sfdry  = 0.215
    mu_sfbet  = 19950.0
    mu_epcap  = 1.321
    mu_fmoil  = 0.823
    mu_floil  = 0.295
    mu_epoil  = 3.827
    mu = np.array([mu_fmmob, mu_sfdry, mu_sfbet, mu_epcap, mu_fmoil, mu_floil, mu_epoil])

    # Coefficient of variation: delta = sigma / mu
    delta_constraint = 1 / np.sqrt(2) 

    delta = []
    for i in range(len(mu)):
        sigma_base = mu[i]*(1.2 - 0.8)/(2*np.sqrt(3))   # base sigma from uniform pdf and +-20%
        delta_base = sigma_base / mu[i]
        delta.append(delta_base)
    print('is delta_base < delta_constraint:', delta < delta_constraint)
    
    sigma = mu * delta_base

    print(mu,'\n',sigma,'\n',delta)

    # Set marginals
    InputOpts = {
        "Marginals": [
            {"Name": 'fmmob',
            "Type": "Gamma",
            "Moments": [mu[0], sigma[0]]
            },
            {"Name": 'sfdry',
            "Type": "Beta",
            "Moments": [mu[1], sigma[1], 0, 1]
            },
            {"Name": 'sfbet',
            "Type": "Gamma",
            "Moments": [mu[2], sigma[2]]
            },
            {"Name": 'epcap',
            "Type": "Gamma",
            "Moments": [mu[3], sigma[3]]
            },
            {"Name": 'fmoil',
            "Type": "Beta",
            "Moments": [mu[4], sigma[4], 0, 1]
            },
            {"Name": 'floil',
            "Type": "Beta",
            "Moments": [mu[5], sigma[5], 0, 1]
            },
            {"Name": 'epoil',
            "Type": "Gamma",
            "Moments": [mu[6], sigma[6]]
            }
        ]
    }

    myInput = uq.createInput(InputOpts)
    uq.print(myInput)

    # # Set sensitivity analysis
    # PCESobol = {
    #     'Type': 'Sensitivity',
    #     'Method': 'Sobol',
    #     'Sobol': {
    #         'Order': 5
    #     }
    # }

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
        PCESobol = {
            'Type': 'Sensitivity',
            'Method': 'Sobol',
            'Sobol': {
                'Order': 3
            }
        }
        QoI_1_SA_CumProd_o.append(uq.createAnalysis(PCESobol))
        # print(QoI_1_SA_CumProd_o)
        print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_1_SA_CumProd_o, PCESobol, idx = -1)}')

        # # QoI 1 - Cumulative oil production
        # print('QoI 1 - Cumulative oil production')
        # myPCE, scaler_X, scaler_Y = eval_PCE(uq, 
        #                                      X = X_ED, 
        #                                      Y = QoI_1_SIM_CumProd_o,
        #                                      training_idx = index, 
        #                                      QoI_col = i
        #                                     )
        # QoI_1_PCE_CumProd_o.append(myPCE)
        # uq.print(QoI_1_PCE_CumProd_o[-1])
        # check_pceAcc(uq,
        #              pce_model = QoI_1_PCE_CumProd_o[-1],
        #              X_SIM = X_ED.iloc[index:, :].to_numpy().tolist(),
        #              Y_SIM = QoI_1_SIM_CumProd_o.iloc[index:,i],
        #              saveDir=experiment_dir + '/PCEsGraphs',
        #              fileName = 'QoI_1_CumProd_o.png',
        #              scaler_X = scaler_X,
        #              scaler_Y = scaler_Y
        # )
        # QoI_1_SA_CumProd_o.append(uq.createAnalysis(PCESobol))
        # print(QoI_1_SA_CumProd_o)
        # print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_1_SA_CumProd_o, PCESobol, idx = -1)}')
        # exit()

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
        PCESobol = {
            'Type': 'Sensitivity',
            'Method': 'Sobol',
            'Sobol': {
                'Order': 3
            }
        }
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
        PCESobol = {
            'Type': 'Sensitivity',
            'Method': 'Sobol',
            'Sobol': {
                'Order': 3
            }
        }
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
        PCESobol = {
            'Type': 'Sensitivity',
            'Method': 'Sobol',
            'Sobol': {
                'Order': 3
            }
        }
        QoI_2_SA_PressDrop.append(uq.createAnalysis(PCESobol))
        print(QoI_2_SA_PressDrop)
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
        PCESobol = {
            'Type': 'Sensitivity',
            'Method': 'Sobol',
            'Sobol': {
                'Order': 3
            }
        }
        QoI_3_SA_OilRecFac.append(uq.createAnalysis(PCESobol))
        print(f'Sum of all {PCESobol['Sobol']['Order']} order Sobol indexes = {SumSobolIndex(QoI_3_SA_OilRecFac, PCESobol, idx = -1)}')

    # Plot Sobol indixes
    nOrder = 3
    plot_SobolIndex(nOrder, QoI_1_SA_CumProd_o, saveDir=experiment_dir + '/SA_analysis/QoI_1_CumProd_o')
    plot_SobolIndex(nOrder, QoI_1_SA_CumProd_w, saveDir=experiment_dir + '/SA_analysis/QoI_1_CumProd_w')
    plot_SobolIndex(nOrder, QoI_1_SA_CumProd_g, saveDir=experiment_dir + '/SA_analysis/QoI_1_CumProd_g')
    plot_SobolIndex(nOrder, QoI_2_SA_PressDrop, saveDir=experiment_dir + '/SA_analysis/QoI_2_SA_PressDrop')
    plot_SobolIndex(nOrder, QoI_3_SA_OilRecFac, saveDir=experiment_dir + '/SA_analysis/QoI_3_SA_OilRecFac')
