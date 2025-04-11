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
import seaborn as sns
import scienceplots
import plotly.express as px

from uqpylab import sessions, display_util

from aux_funcs import *

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams.update({'font.size': 16})

    # Inputs
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap'
    nSamples = 1493
    ti = 0
    tf = 40000
    write_interval = 200
    times2PCE = [39800] 
    tD = 0.36

    # Get samples from OpenFOAM simulations parameters
    experiment_dir = set_experimentDir(experiment_name)
    X_ED = pd.read_csv(experiment_dir + '/X_ED.csv',header=None)
    indexes_to_remove = [49,476,733,953,914,1443,1451]
    X_ED = X_ED.drop(indexes_to_remove).reset_index(drop=True)
    X_ED = X_ED[:nSamples]
    # print(X_ED.shape)

    # Initialize QoIs array from simulations
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    idxs_times2PCE = [np.argmin(np.abs(t-time)) for time in times2PCE]
    # print(idxs_times2PCE)
    QoI_1_SIM_CumProd_o = pd.read_csv(experiment_dir + '/QoI_1_cumulativeProdOil.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_2_SIM_PressDrop = pd.read_csv(experiment_dir + '/QoI_2_pressDrop.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_3_SIM_OilRecFac = pd.read_csv(experiment_dir + '/QoI_3_OilRecFac.csv', header=None).iloc[:nSamples,idxs_times2PCE]
    QoI_5_SIM_OilBank_length = pd.read_csv(experiment_dir + '/OilBank_length_OilSat_'+str(tD)+'.csv', header=None).iloc[:nSamples,:]
    QoI_5_SIM_OilBank_height = pd.read_csv(experiment_dir + '/OilBank_height_OilSat_'+str(tD)+'.csv', header=None).iloc[:nSamples,:]
    QoI_5_SIM_OilBank_area = pd.read_csv(experiment_dir + '/OilBank_area_OilSat_'+str(tD)+'.csv', header=None).iloc[:nSamples,:]


    # Reshape the dataframe
    X_ED['delta'] = 1 - X_ED.iloc[:,5]/X_ED.iloc[:,4]
    X_ED = pd.concat([X_ED, QoI_1_SIM_CumProd_o, QoI_2_SIM_PressDrop, QoI_3_SIM_OilRecFac, QoI_5_SIM_OilBank_length, QoI_5_SIM_OilBank_height, QoI_5_SIM_OilBank_area], axis=1)
    column_names = ['fmmob', 'SF', 'sfbet', 'epcap', 'fmoil', 'floil', 'epoil', 'delta', 'CumProd_o', 'PressDrop', 'OilRecFac', 'OilBank_length', 'OilBank_height', 'OilBank_area']
    X_ED.columns = column_names
    print(X_ED.head())

    QoIs = ['CumProd_o', 'PressDrop', 'OilRecFac', 'OilBank_length', 'OilBank_height', 'OilBank_area']

    # ##############################
    # PLOTS - ORIGINAL DATAFRAME
    # ##############################

    # Correlation matrix
    for i in range(len(QoIs)):
        corr_toPlot = QoIs[i]
        QoIs_new = [QoIs[j] for j in range(len(QoIs)) if j != i]
        X_ED_new = X_ED.drop(columns=QoIs_new)
        corr_matrix = X_ED_new.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", linewidths=0.5, vmin=-1, vmax=1)
        plt.title('Correlation matrix - QoI: ' + corr_toPlot)
        plt.tight_layout()
        plt.savefig(experiment_dir + '/original_CorrMatrix_' + corr_toPlot + '.png', dpi = 200)
  
    # Scatter plot
    plt.rcParams.update({'font.size': 22})
    for i in range(len(QoIs)):
        corr_toPlot = QoIs[i]
        QoIs_new = [QoIs[j] for j in range(len(QoIs)) if j != i]
        X_ED_new = X_ED.drop(columns=QoIs_new)
        plt.figure(figsize=(10, 8))
        sns.pairplot(X_ED_new, corner=True, kind="hist")
        plt.tight_layout()
        plt.savefig(experiment_dir + '/original_Scatter_' + corr_toPlot + '.png', dpi = 200)

        plt.close()
        g = sns.PairGrid(X_ED_new, x_vars=np.concatenate((['fmmob', 'SF', 'sfbet', 'epcap', 'fmoil', 'floil', 'epoil', 'delta'],[corr_toPlot])), y_vars=[corr_toPlot], height=3.0, aspect=1.0)
        g.map_diag(sns.histplot)
        g.map_offdiag(sns.histplot)
        for ax in g.axes[-1, :]: 
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
        plt.tight_layout()
        g.savefig(experiment_dir + '/hist_QoI_only_' + corr_toPlot + '.png', dpi = 200)
        # plt.show()
        
    plt.rcParams.update({'font.size': 16})
    exit()

    # Relation
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    # X_ED.plot.scatter(x='epoil', y='delta', c='OilRecFac', cmap='nipy_spectral', colorbar=True, alpha=0.5, ax=axes[2])
    X_ED.plot.scatter(x='floil', y='fmoil', c='OilRecFac', cmap='nipy_spectral', colorbar=True, alpha=0.5, ax=axes[0])
    X_ED.plot.scatter(x='epcap', y='fmmob', c='OilRecFac', cmap='nipy_spectral', colorbar=True, alpha=0.5, ax=axes[1])
    plt.tight_layout()
    plt.savefig(experiment_dir + '/original_relation' + '.png', dpi = 200)
    plt.close()
    
    # ##############################
    # FILTERING DATA
    # ##############################

    X_ED = X_ED[(X_ED['fmmob'] < 51000) | (X_ED['epcap'] > 1.25)]
    X_ED = X_ED[(X_ED['fmmob'] > 45000) | (X_ED['epcap'] < 1.435)]
    X_ED = X_ED[(X_ED['fmoil'] > 0.756)]

    # ##############################
    # PLOTS - FILTERED DATAFRAME
    # ##############################

    # # Correlation matrix
    # for i in range(len(QoIs)):
    #     corr_toPlot = QoIs[i]
    #     QoIs_new = [QoIs[j] for j in range(len(QoIs)) if j != i]
    #     X_ED_new = X_ED.drop(columns=QoIs_new)
    #     corr_matrix = X_ED_new.corr()
    #     plt.figure(figsize=(10, 8))
    #     sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", linewidths=0.5, vmin=-1, vmax=1)
    #     plt.title('Correlation matrix - QoI: ' + corr_toPlot)
    #     plt.tight_layout()
    #     plt.savefig(experiment_dir + '/filtered_CorrMatrix_' + corr_toPlot + '.png', dpi = 200)
  
    # Scatter plot
    plt.rcParams.update({'font.size': 22})
    for i in range(len(QoIs)):
        corr_toPlot = QoIs[i]
        QoIs_new = [QoIs[j] for j in range(len(QoIs)) if j != i]
        X_ED_new = X_ED.drop(columns=QoIs_new)
        # plt.figure(figsize=(10, 8))
        # sns.pairplot(X_ED_new, corner=True, kind="hist")
        # plt.tight_layout()
        # plt.savefig(experiment_dir + '/filtered_Scatter_' + corr_toPlot + '.png', dpi = 200)

        plt.close()
        g = sns.PairGrid(X_ED_new, x_vars=np.concatenate((['fmmob', 'SF', 'sfbet', 'epcap', 'fmoil', 'floil', 'epoil', 'delta'],[corr_toPlot])), y_vars=[corr_toPlot], height=3.0, aspect=1.0)
        g.map_diag(sns.histplot)
        g.map_offdiag(sns.histplot)
        for ax in g.axes[-1, :]: 
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
        plt.tight_layout()
        g.savefig(experiment_dir + '/hist_QoI_only_' + corr_toPlot + '.png', dpi = 200)
        plt.show()

    plt.rcParams.update({'font.size': 16})

    # Relation
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    # X_ED.plot.scatter(x='epoil', y='delta', c='OilRecFac', cmap='nipy_spectral', colorbar=True, alpha=0.5, ax=axes[2])
    X_ED.plot.scatter(x='floil', y='fmoil', c='OilRecFac', cmap='nipy_spectral', colorbar=True, alpha=0.5, ax=axes[0])
    X_ED.plot.scatter(x='epcap', y='fmmob', c='OilRecFac', cmap='nipy_spectral', colorbar=True, alpha=0.5, ax=axes[1])
    plt.tight_layout()
    plt.savefig(experiment_dir + '/filtered_relation' + '.png', dpi = 200)
    plt.close()
  
