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


if __name__ == '__main__':
    
    plt.style.use('science')
    # plt.rcParams.update({'font.size': 28})

    Swc = 0.0 # 0.197
    Sgr = 0.0 # 0.013
    Sor = 0.0 # 0.103  

    ExeTime_OF = 2.92
    ExeTime_FOSSIL = 23.5

    tD = 1.0
    files_dir = os.getcwd() + '/Tang/'
    post_filename = '_PVI_' + str(tD) + '.csv'


    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.size'] = 20

    ti = 0
    tf = 4400
    write_interval = 25
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    print(t)

    # df = pd.DataFrame([t.astype(int)]) 
    # df.to_csv('output.csv', index=False, header=False)

    if tD == 0.44:  t_find = 1050
    if tD == 1.00:  t_find = 2425
    if tD == 10.40: t_find = 99950



    # sup_data_dir = os.getcwd() + '/sup_data.csv'
    # sup_data = pd.read_csv(sup_data_dir)
    # x = sup_data['Points:2'][:-1].to_numpy()
    # L = x[-1]
    # x = x/L

    
    experiment_name = ''
    experiment_dir = '.'
    sample_dir = experiment_dir + '/' + experiment_name 
    # print(sample_dir)
    data_dict = pd.read_csv('profiles3D'+str(tD)+'.csv')


    # print(data_dict.shape)
    # print(list(data_dict.columns.values))

    Sa = data_dict['Sa']
    Sb = data_dict['Sb']
    Sc = 1.0 - Sa - Sb

    # Sa = S_norm(Sa,Sgr,Swc,Sgr,Sor)
    # Sb = S_norm(Sb,Swc,Swc,Sgr,Sor)
    # Sc = S_norm(Sc,Sor,Swc,Sgr,Sor)
    print(data_dict.head())
    x = data_dict['Points:2'].to_numpy()
    L = x[-1]
    x = x/L

    fig, axes = plt.subplots(1, 2, figsize=(20, 4))
    
    # CT data
    data_So = pd.read_csv(files_dir + 'So' + post_filename)
    x_SoCT = data_So.iloc[:,0].to_numpy()
    So_CT = data_So.iloc[:,1].to_numpy()

    data_Sg = pd.read_csv(files_dir + 'Sg' + post_filename)
    x_SgCT = data_Sg.iloc[:,0].to_numpy()
    Sg_CT = data_Sg.iloc[:,1].to_numpy()

    data_Sw = pd.read_csv(files_dir + 'Sw' + post_filename)
    x_SwCT = data_Sw.iloc[:,0].to_numpy()
    Sw_CT = data_Sw.iloc[:,1].to_numpy()

    if tD != 10.4:
        # DARTS data
        post_filename = '_PVI_' + str(tD) + '_DARTS.csv'

        data_So = pd.read_csv(files_dir + 'So' + post_filename)
        x_SoDARTS = data_So.iloc[:,0].to_numpy()
        So_DARTS = data_So.iloc[:,1].to_numpy()

        data_Sg = pd.read_csv(files_dir + 'Sg' + post_filename)
        x_SgDARTS = data_Sg.iloc[:,0].to_numpy()
        Sg_DARTS = data_Sg.iloc[:,1].to_numpy()

        data_Sw = pd.read_csv(files_dir + 'Sw' + post_filename)
        x_SwDARTS = data_Sw.iloc[:,0].to_numpy()
        Sw_DARTS = data_Sw.iloc[:,1].to_numpy()

        axes[0].plot(x_SoDARTS,So_DARTS,'-r', linewidth = 2,label='DARTS - (Lyu, X., 2021)')
        # axes[1].plot(x_SwDARTS,Sw_DARTS,'-r', linewidth = 2,label='DARTS - (Lyu, X., 2021)')
        axes[1].plot(x_SgDARTS,Sg_DARTS,'-r', linewidth = 2,label='DARTS - (Lyu, X., 2021)')
        
    # axes[0].scatter(x_SoCT,So_CT,c='k',label='CT scan')
    axes[0].plot(x_SoCT,So_CT,'-ko',label='CT scan - (Tang, J., 2019)')
    # axes[0].plot(x_SoDARTS,So_DARTS,'-r',label='DARTS - (Lyu, X., 2021)')
    # axes[0].plot(x_FOSSIL,So_FOSSIL,'-g',label='FOSSIL')
    axes[0].set_xlabel(r'$x/L$ [-]')
    axes[0].set_ylabel(r'$S_o$ [-]')
    # axes[0].legend(loc='lower right')
    axes[0].set_ylim([0,1])
    axes[0].grid()

    # # axes[1].scatter(x_SwCT,Sw_CT,c='k',label='CT scan')
    # axes[1].plot(x_SwCT,Sw_CT,'-ko',label='CT scan - (Tang, J., 2019)')
    # # axes[1].plot(x_SwDARTS,Sw_DARTS,'-r',label='DARTS - (Lyu, X., 2021)')
    # # axes[1].plot(x_FOSSIL,Sw_FOSSIL,'-g',label='FOSSIL')
    # axes[1].set_xlabel(r'$x/L$ [-]')
    # axes[1].set_ylabel(r'$S_w$ [-]')
    # # axes[1].legend()
    # axes[1].set_ylim([0,1])
    # axes[1].grid()

    # axes[2].scatter(x_SgCT,Sg_CT,c='k',label='CT scan')
    axes[1].plot(x_SgCT,Sg_CT,'-ko',label='CT scan - (Tang, J., 2019)')
    # axes[2].plot(x_SgDARTS,Sg_DARTS,'-r',label='DARTS - (Lyu, X., 2021)')
    # axes[2].plot(x_FOSSIL,Sg_FOSSIL,'-g',label='FOSSIL')
    axes[1].set_xlabel(r'$x/L$ [-]')
    axes[1].set_ylabel(r'$S_g$ [-]')
    axes[1].set_ylim([0,1])
    axes[1].grid()

    axes[0].plot(x[:-1],Sc[:-1],c='b', linewidth = 2, label='ImpesFOAM')
    # axes[0].plot(x_FOSSIL,So_FOSSIL,c='g', linewidth = 2, label='FOSSIL')

    # axes[1].plot(x[:-1],Sb[:-1],c='b', linewidth = 2, label='ImpesFOAM')
    # axes[1].plot(x_FOSSIL,Sw_FOSSIL,c='g', linewidth = 2, label='FOSSIL')

    axes[1].plot(x[:-1],Sa[:-1],c='b', linewidth = 2, label='ImpesFOAM')
    # axes[2].plot(x_FOSSIL,Sg_FOSSIL,c='g', linewidth = 2, label='FOSSIL')

    axes[1].legend(loc='upper right',fontsize=16)

    plt.rcParams['font.size'] = 14
    # fig.suptitle(f'PVI:${tD}$' + '\n\nIMPES-OpenFOAM, nCells = '+str(len(x))+r', ExecutionTime$\approx$'+str(ExeTime_OF) + ' min' + '\nFOSSIL, nCells = '+str(int(len(x_FOSSIL)/4))+r', ExecutionTime$\approx$'+str(ExeTime_FOSSIL)+' min')   

    fig.tight_layout()

    # plt.savefig(f'PVI_{tD}.png', dpi=500)
    plt.savefig(f'PVI_{tD}.pdf', dpi=300)
    plt.show()
    plt.close()


    ps = data_dict['p']
    pD = []
    for i in range(len(ps)):
        p = ps.iloc[i]
        pD.append(p[0] - p[-1])
        # plt.plot(x,ps.iloc[i])
        # plt.show()
        # plt.close()
    plt.plot(t[1:],pD)
    plt.show()