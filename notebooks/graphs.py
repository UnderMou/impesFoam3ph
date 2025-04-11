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
from scipy import stats as st

from aux_funcs import *

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams.update({'font.size': 16})

    # Inputs
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_newrange'
    nSamples = 906 # 1493
    ti = 0
    tf = 40000
    write_interval = 200
    nPCEs = 1
    trainingSet_size = 0.75
    times2PCE = [39800] # np.linspace(2000,tf,nPCEs).astype(int)

    # Get samples from OpenFOAM simulations parameters
    experiment_dir = set_experimentDir(experiment_name)
    X_ED = pd.read_csv(experiment_dir + '/X_ED.csv',header=None)
    # indexes_to_remove = [49,476,733,953,914,1443,1451]
    # X_ED = X_ED.drop(indexes_to_remove).reset_index(drop=True)
    X_ED = X_ED[:nSamples]
    print(X_ED.shape)
    print(X_ED.head())
    # X_ED.iloc[:, 5] = -(X_ED.iloc[:,5]/X_ED.iloc[:,4] - 1)
    X_ED['delta'] = -(X_ED.iloc[:,5]/X_ED.iloc[:,4] - 1)
    column_names = ['fmmob', 'SF', 'sfbet', 'epcap', 'fmoil', 'floil', 'epoil', 'delta']
    X_ED.columns = column_names
    print(X_ED.head())
    
    # Pairplot
    plt.rcParams.update({'font.size': 25})
    g = sns.pairplot(X_ED, corner=True, kind="hist")
    for ax in g.axes[-1, :]:  # Última linha de eixos (onde os labels de x aparecem)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
    plt.savefig(experiment_dir + '/pairplot.pdf', dpi=200)
    plt.tight_layout()
    plt.rcParams.update({'font.size': 16})

    index = round(X_ED.shape[0]*trainingSet_size)
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    idxs_times2PCE = [np.argmin(np.abs(t-time)) for time in times2PCE]
    print(idxs_times2PCE)

    # Initialize QoIs array from simulations
    CumProd_o = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdOil.csv', header=None).iloc[:nSamples,:]
    CumProd_w = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdWater.csv', header=None).iloc[:nSamples,:]
    CumProd_g = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdGas.csv', header=None).iloc[:nSamples,:]
    PressDrop = pd.read_csv(experiment_dir + '/QoI_2_pressDrop.csv', header=None).iloc[:nSamples,:]
    OilRecFac = pd.read_csv(experiment_dir + '/QoI_3_OilRecFac.csv', header=None).iloc[:nSamples,:]
    OilCut = pd.read_csv(experiment_dir + '/QoI_4_OilCut.csv', header=None).iloc[:nSamples,:]

    # Plot pairplot
    
    
    # Plot histograms
    idt = -1 # steady state index
    plot_QoIs_histograms_steadyState(idt,
                                     experiment_dir, 
                                     CumProd_o.to_numpy(),
                                     CumProd_w.to_numpy(),
                                     CumProd_g.to_numpy(),
                                     PressDrop.to_numpy(),
                                     OilRecFac.to_numpy(),
                                     nBins = 30)
    plt.close()

    # Filtering samples
    desiredTime = 39800
    QoI = PressDrop
    QoI_name = 'PressDrop'

    id_desired = np.argmin(np.abs(t-desiredTime))
    # last breakthrough
    id_1 = np.argmin(QoI.iloc[:,id_desired])
    # highest production
    id_2 = np.argmax(QoI.iloc[:,id_desired])
    # Mode on desired time
    Median = np.median(QoI.iloc[:,id_desired])
    id_3 = np.argmin(np.abs(Median - QoI.iloc[:,id_desired]))

    plt.figure(figsize=(7,7))
    # Plot fill bewtween with lines
    percentil_5 = np.zeros(QoI.shape[1])
    percentil_95 = np.zeros(QoI.shape[1])
    print(percentil_5.shape)
    for i in range(QoI.shape[1]):
        percentil_5[i] = np.percentile(QoI.iloc[:,i], 5) 
        percentil_95[i] = np.percentile(QoI.iloc[:,i], 95) 
    # plt.fill_between(t[1:], percentil_5, percentil_95, color='cyan', alpha=0.3, label='Intervalo 5%-95%')

    for i in range(QoI.shape[0]):
        plt.plot(t[1:],QoI.iloc[i,:],c='gray',alpha=0.1)
    # plt.plot(t[1:],percentil_5,'k-.o', markersize=3)
    # plt.plot(t[1:],percentil_95,'k-.o', markersize=3)
    plt.plot(t[1:],percentil_5,'k-.',label='5/95 percentiles')
    plt.plot(t[1:],percentil_95,'k-.')

    plt.plot(t[1:],QoI.iloc[0,:],c='gray',alpha=0.1, label='Samples')

    plt.plot(t[1:],QoI.iloc[id_1,:],c='b',label='Minor '+QoI_name+' on desired time')
    plt.plot(t[1:],QoI.iloc[id_2,:],c='r',label='Major '+QoI_name+' on desired time')
    plt.plot(t[1:],QoI.iloc[id_3,:],c='yellow',label='Median of '+QoI_name+' on desired time')

    plt.plot([desiredTime,desiredTime], [0,1.5*np.max(QoI.iloc[:,-1])], 'k--')
    plt.legend()
    plt.grid()
    plt.savefig(experiment_dir + '/'+QoI_name+'_' + str(desiredTime) + '.png', dpi=200)
    # plt.show()
