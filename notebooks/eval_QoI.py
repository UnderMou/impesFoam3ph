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
from sampleMu import *
from uqpylab import sessions, display_util

def generate_color_vector(N, colormap_name='viridis'):
    """
    Generate a vector of colors based on an integer N.

    Parameters:
        N (int): Number of colors to generate.
        colormap_name (str): Name of the matplotlib colormap to use.
        
    Returns:
        list: A list of RGB tuples representing the colors.
    """
    # Get the colormap
    colormap = plt.get_cmap(colormap_name)
    
    # Generate evenly spaced values in [0, 1] for N colors
    color_values = np.linspace(0, 1, N)
    
    # Map the values to RGB colors using the colormap
    colors = [colormap(value) for value in color_values]
    
    return colors

def S_norm(S,Sr,Swc,Sgr,Sor):
    return np.divide(S-Sr,1-Swc-Sor-Sgr)


if __name__ == '__main__':

    Swc = 0.197
    Sgr = 0.013
    Sor = 0.103


    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.size'] = 14

    ti = 0
    tf = 7000
    write_interval = 250
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    print(t)
    t_find = 1500 # tD = 0.36
    idx = np.argmin(np.abs(t-t_find))
    print(idx)

    sup_data_dir = os.getcwd() + '/sup_data.csv'
    sup_data = pd.read_csv(sup_data_dir)
    x = sup_data['Points:2'][:-1].to_numpy()
    L = x[-1]
    x = x/L

    
    experiment_name = 'SA_CTscan_Fdry_baseCase'
    experiment_dir = '/home/anderson/OpenFOAM/anderson-9/run/impesFoam3ph/simulator/experiments'
    experiment_dir += '/' + experiment_name 
    new_dir = "../simulator/experiments/" + experiment_name

    folders = [name for name in os.listdir(experiment_dir) if os.path.isdir(os.path.join(experiment_dir, name))]
    sorted_folders = sorted(folders, key=lambda x: int(x.split('_')[1]))
    folders = sorted_folders
    print(folders)




    X_ED = pd.read_csv(experiment_dir + '/X_ED.csv',header=None)
    print(X_ED.head())

    # #define format for subplots (1 row and 3 columns)
    # fig, axis = plt.subplots(1, 3)

    # #create histogram for each column in DataFrame
    # X_ED.hist(ax=axis)
    # plt.show()
    # plt.close()

    mask = [True] * len(X_ED)

        
    filter = mask
    print(filter)
    print(X_ED[filter])

    filter_idx = np.where(filter)[0].tolist()
    filter_idx = np.array(filter_idx)
    new_folders = np.array(folders)[filter_idx]

    nsamples = len(new_folders)
    time_steps = len(t)-1
    t = t[:-1]
    prod_w = np.zeros((nsamples, time_steps))
    prod_o = np.zeros((nsamples, time_steps))
    prod_g = np.zeros((nsamples, time_steps))
    press_drop = np.zeros((nsamples, time_steps))

    for k, sample in tqdm(enumerate(new_folders), desc='Reading simulation data'):
        sample_dir = experiment_dir + '/' + sample
        data_dict = parse_openfoam_case(sample_dir, variables=['Fa','Fb','p'])
        aux_df = data_dict.copy()
        aux_df.rename(columns={'Fa': 'Fg', 'Fb': 'Fw'}, inplace=True)
        aux_df['Fo'] = 1 - aux_df['Fg'] - aux_df['Fw']
        # print(aux_df['p'].iloc[1])
        # print(data_dict['p'][:])
        # print(data_dict['p'].iloc[k][0])
        ut = 1.61572e-05

        get = -2
        for col in aux_df.columns:
            if col == 'p':
                aux_df[col] = aux_df[col].apply(lambda x: x[0]) - aux_df[col].apply(lambda x: x[get])
            else:
                aux_df[col] = aux_df[col].apply(lambda x: x[get])
        # print(len(aux_df['p']))
        # exit()


        for j in range(1, time_steps):
            # print(j)
            prod_w[k, j] = ut*np.trapz(aux_df['Fw'].iloc[:j], dx=write_interval)
            prod_o[k, j] = ut*np.trapz(aux_df['Fo'].iloc[:j], dx=write_interval)
            prod_g[k, j] = ut*np.trapz(aux_df['Fg'].iloc[:j], dx=write_interval)
            press_drop[k, j] = aux_df['p'].iloc[j]

    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    for i in range(nsamples):

        axes[0,0].plot(t,prod_o[i,:],c='k',label=str(i))

        axes[0,1].plot(t,prod_w[i,:],c='b',label=str(i))

        axes[1,0].plot(t,prod_g[i,:],c='g',label=str(i))
        

    axes[0,0].set_xlabel(r'$t$ [s]')
    axes[0,0].set_ylabel(r'Oil production [PV]')
    # axes[0].legend(loc='lower right')
    axes[0,0].grid()

    axes[0,1].set_xlabel(r'$t$ [s]')
    axes[0,1].set_ylabel(r'Water production [PV]')
    # axes[1].legend()
    axes[0,1].grid()

    axes[1,0].set_xlabel(r'$t$ [s]')
    axes[1,0].set_ylabel(r'Gas production [PV]')
    # axes[2].legend()
    axes[1,0].grid()
    # plt.savefig(new_dir + '/production_curves.png',dpi=300)
    # plt.show()

    for i in range(nsamples):
        axes[1,1].plot(t,press_drop[i,:],c='r')
    axes[1,1].set_xlabel(r'$t$ [s]')
    axes[1,1].set_ylabel(r'$\Delta p$ [Pa]')
    axes[1,1].grid()
    plt.savefig(new_dir + '/pressure_drop.png',dpi=300)
    plt.show()
    # plt.close()

    # ############ UQ[py]Lab - PCE_model #############

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

    # # Reference values
    # sfdryR = 0.215
    # fmmobR = 50000.0
    # sfbetR = 19950.0    

    # dpr = 0.2
    # dprfmmob = 1e-4

    # dists = {
    # 'num_vars': 3,
    # 'names': ['sfdry', 'fmmob', 'sfbet'],
    # 'bounds': [[sfdryR*(1-dpr), sfdryR*(1+dpr)],
    #             [fmmobR*(1-dprfmmob), fmmobR*(1+dprfmmob)],
    #             [sfbetR*(1-dpr), sfbetR*(1+dpr)]]
    # }


    # Reference values
    sfdryR = 0.215
    fmmobR = 50000.0
    sfbetR = 19950.0    

    dpr = 0.2

    dists = {
    'num_vars': 3,
    'names': ['sfdry', 'fmmob', 'sfbet'],
    'bounds': [[sfdryR*(1-dpr), sfdryR*(1+dpr)],
                [fmmobR*(1-dpr), fmmobR*(1+dpr)],
                [sfbetR*(1-dpr), sfbetR*(1+dpr)]]
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
            }
        ]
    }

    PCESobol = {
        'Type': 'Sensitivity',
        'Method': 'Sobol',
        'Sobol': {
            'Order': 1
        }
    }

    myInput = uq.createInput(InputOpts)

    n_pces = prod_o.shape[1]
    PCEs_prod = []
    PCEs_pres = []

    LOO_prod = []
    LOO_pres = []
    ValErr_prod = []
    ValErr_pres = []

    SA_prod = []
    SA_pres = []

    train_size = 0.75
    index = round(X_ED.shape[0]*train_size)

    # n_pces = 10
    io = 5
    for i in range(io,n_pces):
        print("PCE model no." + str(i) + "/" + str(n_pces))
        print("Time: " + str(t[i]))

        # OIL PRODUCTION
        MetaOpts_prod = {
            'Type': 'Metamodel',
            'MetaType': 'PCE'
        }
        MetaOpts_prod['Degree'] = np.arange(3,5).tolist()
        MetaOpts_prod['ExpDesign'] = {
            'X': X_ED.iloc[0:index,:].to_numpy().tolist(),
            'Y': prod_o[0:index,i].tolist()
        }
        MetaOpts_prod['ValidationSet'] = {
        'X': X_ED.iloc[index:, :].to_numpy().tolist(),
        'Y': prod_o[index:,i].tolist()
        }
        myPCE_prod = uq.createModel(MetaOpts_prod)
        uq.print(myPCE_prod)

        PCEs_prod.append(myPCE_prod)
        LOO_prod.append(myPCE_prod['Error']['LOO'])
        ValErr_prod.append(myPCE_prod['Error']['Val'])

        # SA - Sobol indices PCE
        PCESobolAnalysis = uq.createAnalysis(PCESobol)
        SA_prod.append(PCESobolAnalysis)


        # PRESSURE DROP
        MetaOpts_pres = {
            'Type': 'Metamodel',
            'MetaType': 'PCE'
        }
        MetaOpts_pres['Degree'] = np.arange(3,10).tolist()
        MetaOpts_pres['ExpDesign'] = {
            'X': X_ED.iloc[0:index,:].to_numpy().tolist(),
            'Y': press_drop[0:index,i].tolist()
        }
        MetaOpts_pres['ValidationSet'] = {
            'X': X_ED.iloc[index:, :].to_numpy().tolist(),
            'Y': press_drop[index:,i].tolist()
        }
        myPCE_pres = uq.createModel(MetaOpts_pres)
        uq.print(myPCE_pres)

        PCEs_pres.append(myPCE_pres)
        LOO_pres.append(myPCE_pres['Error']['LOO'])
        ValErr_pres.append(myPCE_pres['Error']['Val'])

        # SA - Sobol indices PCE
        PCESobolAnalysis = uq.createAnalysis(PCESobol)
        SA_pres.append(PCESobolAnalysis)


    # PLOT SOBOL INDICES
    plt.rcParams['font.size'] = 15
    fig, axes = plt.subplots(1, 2, figsize=(15, 9))
    # Oil production
    Si_prod = []
    for SA in SA_prod:
        Si_prod.append(SA['Results']['FirstOrder'])
    df = pd.DataFrame(Si_prod, columns = dists['names'])
    for i in range(df.shape[1]):
        axes[0].scatter(t[io:n_pces],df[dists['names'][i]], label='S'+dists['names'][i]+'-PCE', s=15)
    axes[0].set_xlabel(r'$t$ [s]')
    axes[0].set_ylabel(r'$S_i$ Oil production [-]')
    axes[0].grid(True)
    axes[0].legend()
    
    # Pressure drop
    Si_pres = []
    for SA in SA_pres:
        Si_pres.append(SA['Results']['FirstOrder'])
    df = pd.DataFrame(Si_pres, columns = dists['names'])
    for i in range(df.shape[1]):
        axes[1].scatter(t[io:n_pces],df[dists['names'][i]], label='S'+dists['names'][i]+'-PCE', s=15)
    axes[1].set_xlabel(r'$t$ [s]')
    axes[1].set_ylabel(r'$S_i$ Presure drop [-]')
    axes[1].grid(True)
    axes[1].legend()

    plt.savefig(new_dir + '/SobolIndices.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()