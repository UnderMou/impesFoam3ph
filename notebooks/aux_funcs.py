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
import csv
import seaborn as sns
from scipy.stats import qmc
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import scipy.integrate as integrate
import matplotlib.patches as patches

def read_openfoam_field(file_path):
    """
    Reads an OpenFOAM field file and returns the data as a NumPy array.

    Parameters:
        file_path (str): Path to the OpenFOAM field file.

    Returns:
        np.ndarray: NumPy array with the field data.
    """
    try:
        with open(file_path, 'r') as f:
            content = f.readlines()
        
        # Find the 'internalField' line, and following line is the number of elements
        start_index = next(i for i, line in enumerate(content) if line.startswith('internalField'))
        num_elements = int(content[start_index + 1])
        
        # Extract the data block
        data = content[start_index + 3:start_index + 3 + num_elements]
        
        # Parse data into NumPy array
        values = []
        for line in data:
            line = line.strip().strip('()')
            if ' ' in line:  # Vector or multiple values
                try:
                    values.append(np.array([float(x) for x in line.split()]))
                except ValueError:
                    print(f"Warning: Skipping malformed vector line: {line}")
            else:  # Single value
                try:
                    values.append(float(line))
                except ValueError:
                    print(f"Warning: Skipping malformed scalar line: {line}")

        return np.array(values)

    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def parse_openfoam_case(case_dir, variables=['p', 'Sa', 'Sb', 'U', 'Ua', 'Ub', 'Uc', 'Fa', 'Fb'], time_dirs=None):
    """
    Parses the OpenFOAM case directory structure and reads all field data.
        XXX Note: The default list of variables is too expensive for large samples.

    Parameters:
        case_dir (str): Path to the root directory of the OpenFOAM case.
        variables (list): List of field names to read. 
            Default is pressure ('p'), Saturations ('Sa', 'Sb'), 
            total and velocities ('U', 'Ua', 'Ub', 'Uc') and phase fluxes ('Fa', 'Fb').
        
    Returns:
        pd.DataFrame: Pandas DataFrame with the field data, where each column is a variable 
            and each row is a time step. Each cell contains an array with the field data.
    """
    data = {}

    # Iterate over time directories, e.g. '50', '100', '200', ...
    if time_dirs is None:
        time_dirs = sorted([d for d in os.listdir(case_dir) if d.isdigit() and int(d) > 0], key=lambda x: int(x))
    else:
        if type(time_dirs) == str:
            time_dirs = [time_dirs]
        time_dirs = [str(t) for t in time_dirs]

    for time_dir in time_dirs:
        
        time_path = os.path.join(case_dir, time_dir)

        data[time_dir] = {}

        # Iterate over field/var files in the time directory, e.g. 'U', 'p', 'S', ...
        for field_file in variables:
            field_path = os.path.join(time_path, field_file)
            try:
                data[time_dir][field_file] = read_openfoam_field(field_path)
                # print(f"Read {field_file} from {time_dir}")
            except Exception as e:
                print(f"Error reading {field_file} in {time_dir}: {e}")

    # Convert to DataFrame
    data = pd.DataFrame(data)
    data = data.transpose()
    data.index = data.index.astype(int)

    return data

def run_simulation(base_dir, experiment_name, 
                   fmmob=None,
                   sfdry=None,
                   sfbet=None,
                   fmcap=None,
                   epcap=None,
                   fmoil=None,
                   floil=None,
                   epoil=None,
                   verbose=True):
    """
    Runs an OpenFOAM simulation with the given parameters.

    Parameters:
        base_dir (str): Path to the base OpenFOAM case directory.
        experiment_name (str): Name of the experiment to create.
    """
    assert sfdry is not None and fmmob is not None and sfbet is not None, "One or more variables are None"

    new_dir = "../simulator/experiments/" + experiment_name

    try:
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        else:
            if verbose:
                print(" -- The directory already exists. Files will be overwritten. --")
            
        result = subprocess.run(
            ["cp", "-a", f'{base_dir}/.', new_dir],
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error copying the files:", e.stderr)

    env = Environment(
        loader=FileSystemLoader(new_dir),
        trim_blocks=True,
        lstrip_blocks=True
    )

    template = env.get_template('constant/transportProperties') # Location of the template file with foam parameters

    output = template.render(fmmob=fmmob,
                             sfdry=sfdry,
                             sfbet=sfbet,
                             fmcap=fmcap,
                             epcap=epcap,
                             fmoil=fmoil,
                             floil=floil,
                             epoil=epoil
                            )

    # Overwrite the template transportProperties file with the new values
    with open(os.path.join(new_dir, 'constant', 'transportProperties'), 'w') as f:
        f.write(output)

    try:
        result = subprocess.run(
            ["bash", "run_solver.sh", new_dir],  
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error running the script:", e.stderr)

def process_simulation(params, base_case_directory, experiment_name):
    i, (fmmob, sfdry, sfbet, fmcap, epcap, fmoil, floil, epoil) = params
    experiment_name = experiment_name + f"/sample_{i:03d}"
    run_simulation(base_case_directory, experiment_name, 
                    fmmob=fmmob,
                    sfdry=sfdry,
                    sfbet=sfbet,
                    fmcap=fmcap,
                    epcap=epcap,
                    fmoil=fmoil,
                    floil=floil,
                    epoil=epoil, 
                    verbose=False)

def S_norm(S,Sr,Swc,Sgr,Sor):
    return np.divide(S-Sr,1-Swc-Sor-Sgr)

def set_experimentDir(experiment_name):
    experiment_dir  = os.getcwd() + "/../simulator/experiments/" + experiment_name
    return experiment_dir

def get_folders(experiment_dir):
    # folders = [name for name in os.listdir(experiment_dir) if os.path.isdir(os.path.join(experiment_dir, name))]
    folders = [name for name in os.listdir(experiment_dir) if os.path.isdir(os.path.join(experiment_dir, name)) and name != "SA_analysis" and name != "PCEsGraphs"]
    sorted_folders = sorted(folders, key=lambda x: int(x.split('_')[1]))
    folders = sorted_folders
    return folders

def unfold_residuals(residuals):
    return residuals['Swc'], residuals['Sgr'], residuals['Sor']

def get_Saturations(folders, experiment_dir, time_idx, residuals):
    Swc, Sgr, Sor = unfold_residuals(residuals)
    Sas = []
    Sbs = []
    Scs = []
    for i,sample in tqdm(enumerate(folders), desc='Reading simulation data'):
        sample_dir = experiment_dir + '/' + sample
        data_dict = parse_openfoam_case(sample_dir, variables=['Sa','Sb'])

        Sa = data_dict['Sa'].iloc[time_idx]
        Sb = data_dict['Sb'].iloc[time_idx]
        Sc = 1.0 - Sa - Sb

        Sas.append(S_norm(Sa,Sgr,Swc,Sgr,Sor))
        Sbs.append(S_norm(Sb,Swc,Swc,Sgr,Sor))
        Scs.append(S_norm(Sc,Sor,Swc,Sgr,Sor))
    return Sas, Sbs, Scs

def get_samples(experiment_dir):
    X_ED = pd.read_csv(experiment_dir + '/X_ED.csv',header=None)
    return X_ED.to_numpy()

def get_CTscanData(So_CT_file, Sg_CT_file, Sw_CT_file):
    data_So = pd.read_csv(So_CT_file)
    x_SoCT = data_So.iloc[:,0].to_numpy()
    So_CT = data_So.iloc[:,1].to_numpy()

    data_Sg = pd.read_csv(Sg_CT_file)
    x_SgCT = data_Sg.iloc[:,0].to_numpy()
    Sg_CT = data_Sg.iloc[:,1].to_numpy()

    data_Sw = pd.read_csv(Sw_CT_file)
    x_SwCT = data_Sw.iloc[:,0].to_numpy()
    Sw_CT = data_Sw.iloc[:,1].to_numpy()

    exp_data = {
        'x_SoCT': x_SoCT,
        'So_CT': So_CT,
        'x_SgCT': x_SgCT,
        'Sg_CT': Sg_CT,
        'x_SwCT': x_SwCT,
        'Sw_CT': Sw_CT
    }

    return exp_data

def get_DARTSdata(Sg_DARTS_file, Sw_DARTS_file, So_DARTS_file):

    data_Sg_DARTS = pd.read_csv(Sg_DARTS_file)
    x_Sg_DARTS = data_Sg_DARTS.iloc[:,0].to_numpy()
    Sg_DARTS = data_Sg_DARTS.iloc[:,1].to_numpy()

    data_Sw_DARTS = pd.read_csv(Sw_DARTS_file)
    x_Sw_DARTS = data_Sw_DARTS.iloc[:,0].to_numpy()
    Sw_DARTS = data_Sw_DARTS.iloc[:,1].to_numpy()

    data_So_DARTS = pd.read_csv(So_DARTS_file)
    x_So_DARTS = data_So_DARTS.iloc[:,0].to_numpy()
    So_DARTS = data_So_DARTS.iloc[:,1].to_numpy()

    DARTS_data = {
        'x_So_DARTS': x_So_DARTS,
        'So_DARTS': So_DARTS,
        'x_Sg_DARTS': x_Sg_DARTS,
        'Sg_DARTS': Sg_DARTS,
        'x_Sw_DARTS': x_Sw_DARTS,
        'Sw_DARTS': Sw_DARTS
    }

    return DARTS_data

def eval_production(experiment_dir, folders, reservoir, nsamples,time_steps,nAvoid,residuals,write_interval):
    prod_w = np.zeros((nsamples, time_steps))
    prod_o = np.zeros((nsamples, time_steps))
    prod_g = np.zeros((nsamples, time_steps))

    for i, sample in tqdm(enumerate(folders), desc='Evaluation productions from simulation data'):
        sample_dir = experiment_dir + '/' + sample
        data_dict = parse_openfoam_case(sample_dir, variables=['Sa', 'Sb', 'Ua', 'Ub', 'Uc'])
        aux_df = data_dict.copy()
        aux_df['Uc'] = aux_df['Uc'].iloc[:].apply(lambda x : x[nAvoid][2])
        aux_df['Ub'] = aux_df['Ub'].iloc[:].apply(lambda x : x[nAvoid][2])
        aux_df['Ua'] = aux_df['Ua'].iloc[:].apply(lambda x : x[nAvoid][2])
        # aux_df['Sa'] = aux_df['Sa'].iloc[:].apply(lambda x : x[nAvoid])
        # aux_df['Sb'] = aux_df['Sb'].iloc[:].apply(lambda x : x[nAvoid])
        # aux_df['Sc'] = 1.0 - aux_df['Sa'] - aux_df['Sb']
        
        for j in range(1, time_steps):
            prod_w[i, j] = aux_df['Ub'].iloc[j]*write_interval*reservoir['A']
            prod_o[i, j] = aux_df['Uc'].iloc[j]*write_interval*reservoir['A']
            prod_g[i, j] = aux_df['Ua'].iloc[j]*write_interval*reservoir['A']

    # Swc, Sgr, Sor = unfold_residuals(residuals)
    prod_w /= (reservoir['PoreVol'])
    prod_o /= (reservoir['PoreVol'])
    prod_g /= (reservoir['PoreVol'])

    return prod_w, prod_o, prod_g

def eval_cumulativeProduction(prod_w, prod_o, prod_g,nsamples, time_steps):
    cuml_prod_w = np.zeros((nsamples, time_steps))
    cuml_prod_o = np.zeros((nsamples, time_steps))
    cuml_prod_g = np.zeros((nsamples, time_steps))
    for i in tqdm(range(nsamples), desc='Evaluation cumulative production'):
        for j in range(1, time_steps):
            cuml_prod_w[i, j] = np.sum(prod_w[i, :j])
            cuml_prod_o[i, j] = np.sum(prod_o[i, :j])
            cuml_prod_g[i, j] = np.sum(prod_g[i, :j])
    
    return cuml_prod_w, cuml_prod_o, cuml_prod_g

def eval_pressureDrop(experiment_dir, folders, nsamples,time_steps,nAvoid):
    pressDrop = np.zeros((nsamples, time_steps))
    for i, sample in tqdm(enumerate(folders), desc='Evaluation pressure drop from simulation data'):
        sample_dir = experiment_dir + '/' + sample
        data_dict = parse_openfoam_case(sample_dir, variables=['p'])
        for j in range(1, time_steps):
            # Get pressure at outlet an d inlet at time j
            p_out = data_dict['p'].iloc[j][:nAvoid][-1]
            p_in = data_dict['p'].iloc[j][:nAvoid][0]

            pressDrop[i, j] = p_in - p_out

    return pressDrop

def plot_production(nsamples, prod_o, prod_w, prod_g, t, experiment_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i in range(nsamples):
        axes[0].plot(t[:-1],prod_o[i,:],c='k')
        axes[1].plot(t[:-1],prod_w[i,:],c='b')
        axes[2].plot(t[:-1],prod_g[i,:],c='g')
    axes[0].grid()
    axes[0].set_xlabel('Time [s]')
    axes[0].set_ylabel('Oil production [PV/s]')
    axes[1].grid()
    axes[1].set_xlabel('Time [s]')
    axes[1].set_ylabel('Water production [PV/s]')
    axes[2].grid()
    axes[2].set_xlabel('Time [s]')
    axes[2].set_ylabel('Gas production [PV/s]')
    fig.suptitle('Phase production rate')
    plt.tight_layout()
    plt.savefig(experiment_dir + '/QoI_1_phaseProductionRate.png', dpi=500)
    # plt.show() 

def plot_Cumulativeproduction(nsamples, cuml_prod_o, cuml_prod_w, cuml_prod_g, t, experiment_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i in range(nsamples):
        axes[0].plot(t[:-1],cuml_prod_o[i,:],c='k')
        axes[1].plot(t[:-1],cuml_prod_w[i,:],c='b')
        axes[2].plot(t[:-1],cuml_prod_g[i,:],c='g')
    axes[0].grid()
    axes[0].set_xlabel('Time [s]')
    axes[0].set_ylabel('Oil cumulative production [PV]')
    axes[1].grid()
    axes[1].set_xlabel('Time [s]')
    axes[1].set_ylabel('Water cumulative production [PV]')
    axes[2].grid()
    axes[2].set_xlabel('Time [s]')
    axes[2].set_ylabel('Gas cumulative production [PV]')
    fig.suptitle('Phase cumulative production')
    plt.tight_layout()
    plt.savefig(experiment_dir + '/QoI_1_phaseCumulativeProduction.png', dpi=500)
    # plt.show() 

def plot_pressDrop(nsamples, pressDrop, t, experiment_dir):
    fig, axes = plt.subplots(1, 1, figsize=(5, 5))
    for i in range(nsamples):
        axes.plot(t[:-1],pressDrop[i,:],c='r')
    axes.grid()
    axes.set_xlabel('Time [s]')
    axes.set_ylabel('Pressure drop [Pa]')
    plt.tight_layout()
    plt.savefig(experiment_dir + '/QoI_2_pressureDrop.png', dpi=500)
    # plt.show() 

def plot_OilRecFactor(nsamples, ORF, t, experiment_dir):
    fig, axes = plt.subplots(1, 1, figsize=(5, 5))
    for i in range(nsamples):
        axes.plot(t[:-1],ORF[i,:],c='k')
    axes.grid()
    axes.set_xlabel('Time [s]')
    axes.set_ylabel('Oil Recovery Factor [%]')
    plt.tight_layout()
    plt.savefig(experiment_dir + '/QoI_3_OilRecFac.png', dpi=500)
    # plt.show() 

def plot_OilCut(nsamples, Oil_cut, t, experiment_dir):
    fig, axes = plt.subplots(1, 1, figsize=(5, 5))
    for i in range(nsamples):
        axes.plot(t[:-1],Oil_cut[i,:],c='k')
    axes.grid()
    axes.set_xlabel('Time [s]')
    axes.set_ylabel('Oil cut [%]')
    plt.tight_layout()
    plt.savefig(experiment_dir + '/QoI_4_OilCut.png', dpi=500)
    # plt.show() 

def plot_PD_ORF_OC(nsamples, t, experiment_dir, pressDrop, ORF, Oil_cut):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i in range(nsamples):
        axes[0].plot(t[:-1],pressDrop[i,:],c='r')
        axes[1].plot(t[:-1],ORF[i,:],c='k')
        axes[2].plot(t[:-1],Oil_cut[i,:],c='k')
    
    axes[0].grid()
    axes[0].set_xlabel('Time [s]')
    axes[0].set_ylabel('Pressure drop [Pa]')
    axes[1].grid()
    axes[1].set_xlabel('Time [s]')
    axes[1].set_ylabel('Oil Recovery Factor [%]')
    axes[2].grid()
    axes[2].set_xlabel('Time [s]')
    axes[2].set_ylabel('Oil cut [%]')
    plt.tight_layout()
    plt.savefig(experiment_dir + '/QoI_2_3_4.png', dpi=500)

def eval_PCE(uq, X, Y, training_idx, QoI_col):

    MetaOpts = {
            'Type': 'Metamodel',
            'MetaType': 'PCE',
            'Method': 'OMP', # 'OMP', # 'LARS'
            'Degree': np.arange(0,4,1).tolist(),
            # 'SparseDegree': True,
            'ExpDesign': {
                'X': X.iloc[0:training_idx,:].to_numpy().tolist(),
                'Y': Y.iloc[0:training_idx,QoI_col].to_numpy().tolist()
            },
            'ValidationSet': {
                'X': X.iloc[training_idx:, :].to_numpy().tolist(),
                'Y': Y.iloc[training_idx:,QoI_col].to_numpy().tolist()
            }
        }
    myPCE = uq.createModel(MetaOpts)
    return myPCE

def check_pceAcc(uq, pce_model, X_SIM, Y_SIM, saveDir, fileName):
    currentDir = os.getcwd()
    print('Saving PCEs vs. Simulation plot in: ' + saveDir)
    os.makedirs(saveDir, exist_ok=True)
    os.chdir(saveDir)

    plt.figure(figsize=(4,4))
    Y_PCE = uq.evalModel(pce_model, X_SIM)
    plt.scatter(Y_PCE,Y_SIM)
    ypces = [np.min(Y_PCE), np.max(Y_PCE)]
    ysims = [np.min(Y_SIM), np.max(Y_SIM)]
    print("R^2: ", r2_score(Y_SIM, Y_PCE))
    plt.plot(ysims,ysims)
    # plt.grid()
    plt.title(fileName)
    plt.xlabel(r'$QoI_{PCE}$')
    plt.ylabel(r'$QoI_{simulation}$')
    plt.tight_layout()
    # plt.show()
    plt.savefig(fileName, dpi=500)
    # plt.close()
    os.chdir(currentDir)

def plot_QoIs_histograms_steadyState(idt,
                                     experiment_dir, 
                                     cuml_prod_o,
                                     cuml_prod_w,
                                     cuml_prod_g,
                                     pressDrop,
                                     ORF,
                                     nBins = 10):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    axes[0,0].hist(cuml_prod_o[:,idt],bins=nBins)
    # axes[0,0].set_xticklabels([f"{x:.3f}" for x in cuml_prod_o[:, idt]], rotation=45)
    axes[0,0].grid()
    # axes[0,0].set_title('Cumulative oil production [PV]')
    axes[0,0].set_xlabel('Cumulative oil production [PV]')
    
    axes[0,1].hist(cuml_prod_w[:,idt],bins=nBins)
    # axes[0,1].set_xticklabels([f"{x:.3f}" for x in cuml_prod_w[:, idt]], rotation=45)
    axes[0,1].grid()
    # axes[0,1].set_title('Cumulative water production [PV]')
    axes[0,1].set_xlabel('Cumulative water production [PV]')

    axes[0,2].hist(cuml_prod_g[:,idt],bins=nBins)
    # axes[0,2].set_xticklabels([f"{x:.3f}" for x in cuml_prod_g[:, idt]], rotation=45)
    axes[0,2].grid()
    # axes[0,2].set_title('Cumulative gas production [PV]')
    axes[0,2].set_xlabel('Cumulative gas production [PV]')

    axes[1,0].hist(pressDrop[:,idt]/1e6,bins=nBins)
    # axes[1,0].set_xticklabels([f"{x:.3f}" for x in pressDrop[:, idt]/1e6], rotation=45)
    axes[1,0].grid()
    # axes[1,0].set_title('Pressure drop [MPa]')
    axes[1,0].set_xlabel('Pressure drop [MPa]')

    axes[1,1].hist(ORF[:,idt],bins=nBins)
    # axes[1,1].set_xticklabels([f"{x:.3f}" for x in ORF[:, idt]], rotation=45)
    axes[1,1].grid()
    # axes[1,1].set_title('Oil recovery factor [%]')
    axes[1,1].set_xlabel('Oil recovery factor [%]')

    plt.savefig(experiment_dir + '/QoIs_histogram_steady.png', dpi=500)
    plt.suptitle('QoIs histogram at steady state')
    plt.tight_layout()
    # plt.show()

def SumSobolIndex(SA, PCESobol, idx):
    cumSum = 0
    for i in range(PCESobol['Sobol']['Order']):
        cumSum += np.sum(SA[idx]['Results']['AllOrders'][i])
    return cumSum

def generate_combinations(n, r): 
    if r == 1: ans = [str(i+1) for i in range(n)]
    else:
        elements = list(range(1, n + 1))
        ans = [str(combo) for combo in itertools.combinations(elements, r)]  
    return ans

def plot_SobolIndex(nOrder, QoI_SA, saveDir):
    
    print('Saving Sobol indexes in: ' + saveDir)

    os.makedirs(saveDir, exist_ok=True)

    version = 3

    if version == 1:
        for i in range(nOrder):

            n = len(QoI_SA[-1]['Results']['VariableNames'])
            r = i+1

            title = [str(idx+1)+':'+var+' ' for idx,var in enumerate(QoI_SA[-1]['Results']['VariableNames'])]
            
            formatted_labels = generate_combinations(n, r)

            plt.figure(figsize=(10,3))
            plt.bar(formatted_labels,
                    QoI_SA[-1]['Results']['AllOrders'][i])
            plt.xticks(rotation=90)
            plt.ylabel(f'{i+1}-order Sobol index')
            plt.title(title)
            plt.tight_layout()
            plt.savefig(saveDir + f'/{i+1}-order-SobolIndex.png', dpi=500)
            plt.close()

    if version == 2:

        plt.rcParams.update({'font.size': 22})
        fig, axes = plt.subplots(nOrder+1, 1, figsize=(15, 10))
        
        axes[-1].bar(QoI_SA[-1]['Results']['VariableNames'],
                QoI_SA[-1]['Results']['Total'])
        axes[-1].set_xticklabels(QoI_SA[-1]['Results']['VariableNames'], rotation=0)
        axes[-1].set_ylabel(r'$S_T$')
        axes[-1].set_title('Total Sobol index')
        axes[-1].set_ylim([0,1])

        for i in range(nOrder):

            n = len(QoI_SA[-1]['Results']['VariableNames'])
            r = i+1

            title = [str(idx+1)+':'+var+' ' for idx,var in enumerate(QoI_SA[-1]['Results']['VariableNames'])]
            
            formatted_labels = generate_combinations(n, r)

            axes[i].bar(formatted_labels,
                    QoI_SA[-1]['Results']['AllOrders'][i])
            axes[i].set_xticklabels(formatted_labels, rotation=90)
            axes[i].set_ylabel(f'{i+1}-order')
            axes[i].set_title(title)

            axes[i].set_ylim([0,1])
        
        plt.tight_layout()
        plt.savefig(saveDir + f'_SobolIndexes.png', dpi=500)
        plt.close()

    if version == 3:

        plt.rcParams.update({'font.size': 22})
        # fig, axes = plt.subplots(1, 1, figsize=(10, 5))
        plt.figure(figsize=(7,7))
        plt.bar(QoI_SA[-1]['Results']['VariableNames'],
                QoI_SA[-1]['Results']['Total'])
        plt.xticklabels(QoI_SA[-1]['Results']['VariableNames'], rotation=0)
        plt.ylabel(r'$S_T$')
        plt.title('Total Sobol index')
        plt.ylim([0,1])

def remove_indexes(lst, indexes):
    return [value for i, value in enumerate(lst) if i not in set(indexes)]

def eval_oilBank_PulseDuration_fromOilCut(Oil_cut, t):
    t = t[1:]
    OilBank_PulseDur = np.zeros(Oil_cut.shape[0])
    for i in range(Oil_cut.shape[0]):
        OilCut_grad = np.gradient(Oil_cut.iloc[i,:].to_numpy())
        OilCut_grad = np.nan_to_num(OilCut_grad, nan=0.0)
        init_idt = np.argmax(OilCut_grad)
        # print(init_idt, t[init_idt])
        end_idt = np.argmin(OilCut_grad)
        # print(end_idt, t[end_idt])
        OilBank_PulseDur[i] = t[end_idt] - t[init_idt]
        # print(OilBank_PulseDur[i])
    
    return OilBank_PulseDur

def eval_oilBank_PulseHeight_fromOilCut(Oil_cut, t):
    t = t[1:]
    OilBank_PulseHei = np.zeros(Oil_cut.shape[0])
    for i in range(Oil_cut.shape[0]):
        OilBank_PulseHei[i] = np.max(Oil_cut.iloc[i,:])

    return OilBank_PulseHei

def eval_oilBank_AreaBelow_fromOilCut(Oil_cut, t, write_interval):
    t = t[1:-1]
    OilBank_Area = np.zeros(Oil_cut.shape[0])
    for i in range(Oil_cut.shape[0]):
        OilBank_Area[i] = integrate.trapz(Oil_cut.iloc[i,1:].to_numpy() / write_interval, t)
    
    return OilBank_Area

def plot_oilBank(OilBank_PulseDur, OilBank_PulseHei, OilBank_Area):
    fig, axs = plt.subplots(1,3, figsize=(10,5))

    axs[0].hist(OilBank_PulseDur, bins=7)
    axs[0].set_xlabel('Oil bank duration [s]')
    axs[0].set_ylabel('Frequency')
    axs[0].grid()

    axs[1].hist(OilBank_PulseHei, bins=7)
    axs[1].set_xlabel('Oil bank Height [PV/s]')
    axs[1].set_ylabel('Frequency')
    axs[1].grid()

    axs[2].hist(OilBank_Area, bins=10)
    axs[2].set_xlabel('Oil bank Area [PV]')
    axs[2].set_ylabel('Frequency')
    axs[2].grid()

    plt.tight_layout()
    plt.show()

# def eval_oilBank_length_fromOilSat(Scs, x):

#     OilBank_length = np.zeros(Scs.shape[0])
#     for i in range(Scs.shape[0]):
#         OilSat_grad = np.gradient(Scs.iloc[i,:].to_numpy())
#         OilSat_grad = np.nan_to_num(OilSat_grad, nan=0.0)
#         init_idx = np.argmax(OilSat_grad)
#         # print(init_idt, t[init_idt])
#         end_idx = np.argmin(OilSat_grad)
#         # print(end_idt, t[end_idt])
#         OilBank_length[i] = x[end_idx] - x[init_idx]
    
#     return OilBank_length

def eval_oilBank_length_fromOilSat(Scs, x, Soi=0.46, nAvoid=None):

    OilBank_length = np.zeros(Scs.shape[0])
    init_idxs = np.zeros(Scs.shape[0])
    end_idxs = np.zeros(Scs.shape[0])
    for i in range(Scs.shape[0]):
        OilSat_grad = np.gradient(Scs.iloc[i,:].to_numpy())
        OilSat_grad = np.nan_to_num(OilSat_grad, nan=0.0)
        init_idx = 0
        for j in range(len(Scs.iloc[i,:])):
            if Scs.iloc[i,j] >= Soi: 
                init_idx = j
                break

        end_idx = init_idx + 1
        tol = 1e-3
        find = False
        for j in np.arange(end_idx, len(Scs.iloc[i,:]),1):
            if abs(Scs.iloc[i,j] - Soi) <= tol: 
                end_idx = j
                find = True
                break
        if find == False: end_idx = len(Scs.iloc[i,:]) + nAvoid
        
        # Interpolation
        xSubSet = np.linspace(x[init_idx-1],x[init_idx],100)
        SoSubSet = np.linspace(Scs.iloc[i,init_idx-1],Scs.iloc[i,init_idx],100)
        xExact_idx = np.argmin(np.abs(SoSubSet - Soi))
        xExact = xSubSet[xExact_idx]

        # print(xExact, x[init_idx])
        # plt.figure(figsize=(7,7))
        # plt.scatter(xExact,Soi, c='b')
        # plt.scatter(x,Scs.iloc[i,:],s=2)
        # plt.scatter(x[init_idx], Scs.iloc[i,init_idx], c='k')  
        # plt.scatter(x[end_idx], Scs.iloc[i,end_idx], c='k')       

        # plt.scatter(xSubSet, SoSubSet, s=0.1,c='k')
        # plt.plot([xExact,x[end_idx]], [Soi, Soi], 'k--') 
        # plt.show()

        # OilBank_length[i] = x[end_idx] - x[init_idx]
        OilBank_length[i] = x[end_idx] - xExact
        init_idxs[i] = init_idx
        end_idxs[i] = end_idx
    
    return OilBank_length, init_idxs.astype(int), end_idxs.astype(int)

def eval_oilBank_height_fromOilSat(Scs, x, init_idxs, end_idxs, Soi = 0.46):

    OilBank_height = np.zeros(Scs.shape[0])
    for i in range(Scs.shape[0]):

        So_Bank = Scs.iloc[i, init_idxs[i]:end_idxs[i]]

        OilBankSat_mean = np.median(So_Bank)
        OilBank_height[i] = OilBankSat_mean - Soi
    
    return OilBank_height

def eval_oilBank_AreaBelow_fromOilCut(Scs, x, init_idxs, end_idxs, Soi = 0.46):

    OilBank_Area = np.zeros(Scs.shape[0])

    for i in range(Scs.shape[0]):

        x_Bank = x[init_idxs[i]:end_idxs[i]] 
        So_Bank = Scs.iloc[i, init_idxs[i]:end_idxs[i]]
        So_init = Soi * np.ones_like(So_Bank)

        area_up = integrate.trapz(So_Bank, x_Bank)
        area_down = integrate.trapz(So_init, x_Bank)
        OilBank_Area[i] = area_up - area_down
        # print(area_down, area_up, OilBank_Area[i])
    
    return OilBank_Area

def plot_fromOilSat(Scs, x, init_idxs, end_idxs, Soi = 0.46):

    for i in range(Scs.shape[0]):

        x_Bank = x[init_idxs[i]:end_idxs[i]] 
        So_Bank = Scs.iloc[i, init_idxs[i]:end_idxs[i]]
        So_init = Soi * np.ones_like(So_Bank)

        OilBankSat_mean = np.median(So_Bank)

        print(Scs.iloc[i,init_idxs[i]], Scs.iloc[i,end_idxs[i]],OilBankSat_mean, OilBankSat_mean-Soi)
        fig, ax = plt.subplots(figsize=(7,7))
        plt.scatter(x,Scs.iloc[i,:],s=2)
        plt.plot(x,Scs.iloc[i,:],'.-')
        plt.plot(x_Bank,So_Bank,c='k')
        plt.plot(x_Bank,So_init,c='k')

        plt.xlabel(r'$x/L$ [-]')
        plt.ylabel(r'$S_o$ [-]')
        plt.grid()

        # plt.plot([x_Bank[0], x_Bank[-1]],[OilBankSat_mean, OilBankSat_mean],'r--', label='Oil bank height')

        plt.scatter(x[init_idxs[i]], Scs.iloc[i,init_idxs[i]], c='k', label=r'Oil-bank initial/end point')
        plt.scatter(x[end_idxs[i]], Scs.iloc[i,end_idxs[i]], c='k')
        plt.fill_between(x_Bank, So_init, So_Bank, color='b', alpha=0.3, label ='Oil-bank area')
        plt.ylim([0,1])
        plt.legend()
        plt.show()