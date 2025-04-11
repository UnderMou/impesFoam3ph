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

    ExeTime_OF = 2.92
    ExeTime_FOSSIL = 23.5

    tD = 10.4
    files_dir = os.getcwd() + '/Tang/'
    post_filename = '_PVI_' + str(tD) + '.csv'


    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.size'] = 16

    ti = 0
    tf = 40000
    write_interval = 25
    t = np.linspace(ti,tf,int(tf/write_interval + 1))

    # df = pd.DataFrame([t.astype(int)]) 
    # df.to_csv('output.csv', index=False, header=False)

    if tD == 0.44:  t_find = 750
    if tD == 1.00:  t_find = 1705
    if tD == 10.40: t_find = 17725

    # Reading FOSSIL data
    if tD == 0.44:
        FOSSIL_data_PVI044 = pd.read_csv('FOSSIL_CT_PV0.44.csv')
        Sg_FOSSIL_PVI044 = FOSSIL_data_PVI044['Gas_Saturation']
        Sg_FOSSIL = S_norm(Sg_FOSSIL_PVI044,Sgr,Swc,Sgr,Sor)
        So_FOSSIL_PVI044 = FOSSIL_data_PVI044['Oil_Saturation']
        So_FOSSIL = S_norm(So_FOSSIL_PVI044,Sor,Swc,Sgr,Sor)
        Sw_FOSSIL_PVI044 = FOSSIL_data_PVI044['Water_Saturation']
        Sw_FOSSIL = S_norm(Sw_FOSSIL_PVI044,Swc,Swc,Sgr,Sor)
        x_FOSSIL = FOSSIL_data_PVI044['Points:0']

    if tD == 1.00:
        FOSSIL_data_PVI1 = pd.read_csv('FOSSIL_CT_PV1.0.csv')
        Sg_FOSSIL_PVI1 = FOSSIL_data_PVI1['Gas_Saturation']
        Sg_FOSSIL = S_norm(Sg_FOSSIL_PVI1,Sgr,Swc,Sgr,Sor)
        So_FOSSIL_PVI1 = FOSSIL_data_PVI1['Oil_Saturation']
        So_FOSSIL = S_norm(So_FOSSIL_PVI1,Sor,Swc,Sgr,Sor)
        Sw_FOSSIL_PVI1 = FOSSIL_data_PVI1['Water_Saturation']
        Sw_FOSSIL = S_norm(Sw_FOSSIL_PVI1,Swc,Swc,Sgr,Sor)
        x_FOSSIL = FOSSIL_data_PVI1['Points:0']

    if tD == 10.4:
        FOSSIL_data_PVI104 = pd.read_csv('FOSSIL_CT_PV10.4.csv')
        Sg_FOSSIL_PVI104 = FOSSIL_data_PVI104['Gas_Saturation']
        Sg_FOSSIL = S_norm(Sg_FOSSIL_PVI104,Sgr,Swc,Sgr,Sor)
        So_FOSSIL_PVI104 = FOSSIL_data_PVI104['Oil_Saturation']
        So_FOSSIL = S_norm(So_FOSSIL_PVI104,Sor,Swc,Sgr,Sor)
        Sw_FOSSIL_PVI104 = FOSSIL_data_PVI104['Water_Saturation']
        Sw_FOSSIL = S_norm(Sw_FOSSIL_PVI104,Swc,Swc,Sgr,Sor)
        x_FOSSIL = FOSSIL_data_PVI104['Points:0']
    
    L = x_FOSSIL.to_numpy()[-1]
    x_FOSSIL = x_FOSSIL/L
    
    idx = np.argmin(np.abs(t-t_find))
    print(idx)

    # sup_data_dir = os.getcwd() + '/sup_data.csv'
    # sup_data = pd.read_csv(sup_data_dir)
    # x = sup_data['Points:2'][:-1].to_numpy()
    # L = x[-1]
    # x = x/L

    
    experiment_name = ''
    experiment_dir = '/home/anderson/OpenFOAM/anderson-9/run/impesFoam3ph_v1/Lyu_CT_exp'
    sample_dir = experiment_dir + '/' + experiment_name 
    # print(sample_dir)
    data_dict = parse_openfoam_case(sample_dir, variables=['Sa','Sb'])

    # print(data_dict.shape)
    # print(list(data_dict.columns.values))

    Sa = data_dict['Sa'].iloc[idx]
    Sb = data_dict['Sb'].iloc[idx]
    Sc = 1.0 - Sa - Sb

    Sa = S_norm(Sa,Sgr,Swc,Sgr,Sor)
    Sb = S_norm(Sb,Swc,Swc,Sgr,Sor)
    Sc = S_norm(Sc,Sor,Swc,Sgr,Sor)

    x = np.linspace(0,1,len(Sa))
    L = x[-1]
    x = x/L

    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    

    data_So = pd.read_csv(files_dir + 'So' + post_filename)
    x_SoCT = data_So.iloc[:,0].to_numpy()
    So_CT = data_So.iloc[:,1].to_numpy()

    data_Sg = pd.read_csv(files_dir + 'Sg' + post_filename)
    x_SgCT = data_Sg.iloc[:,0].to_numpy()
    Sg_CT = data_Sg.iloc[:,1].to_numpy()

    data_Sw = pd.read_csv(files_dir + 'Sw' + post_filename)
    x_SwCT = data_Sw.iloc[:,0].to_numpy()
    Sw_CT = data_Sw.iloc[:,1].to_numpy()
        
    # axes[0].scatter(x_SoCT,So_CT,c='k',label='CT scan')
    axes[0].plot(x_SoCT,So_CT,'-ko',label='CT scan - (Tang, J., 2019)')
    axes[0].set_xlabel(r'$x/L$ [-]')
    axes[0].set_ylabel(r'$S_o$ [-]')
    # axes[0].legend(loc='lower right')
    axes[0].set_ylim([0,1])
    axes[0].grid()

    # axes[1].scatter(x_SwCT,Sw_CT,c='k',label='CT scan')
    axes[1].plot(x_SwCT,Sw_CT,'-ko',label='CT scan - (Tang, J., 2019)')
    axes[1].set_xlabel(r'$x/L$ [-]')
    axes[1].set_ylabel(r'$S_w$ [-]')
    # axes[1].legend()
    axes[1].set_ylim([0,1])
    axes[1].grid()

    # axes[2].scatter(x_SgCT,Sg_CT,c='k',label='CT scan')
    axes[2].plot(x_SgCT,Sg_CT,'-ko',label='CT scan - (Tang, J., 2019)')
    axes[2].set_xlabel(r'$x/L$ [-]')
    axes[2].set_ylabel(r'$S_g$ [-]')
    axes[2].set_ylim([0,1])
    axes[2].grid()

    axes[0].plot(x[:-1],Sc[:-1],c='b', linewidth = 2, label='IMPES-OpenFOAM')
    axes[0].plot(x_FOSSIL,So_FOSSIL,c='g', linewidth = 2, label='FOSSIL')

    axes[1].plot(x[:-1],Sb[:-1],c='b', linewidth = 2, label='IMPES-OpenFOAM')
    axes[1].plot(x_FOSSIL,Sw_FOSSIL,c='g', linewidth = 2, label='FOSSIL')

    axes[2].plot(x[:-1],Sa[:-1],c='b', linewidth = 2, label='IMPES-OpenFOAM')
    axes[2].plot(x_FOSSIL,Sg_FOSSIL,c='g', linewidth = 2, label='FOSSIL')

    axes[1].legend(loc='upper center')

    plt.rcParams['font.size'] = 14
    fig.suptitle(f'PVI:${tD}$' + '\n\nIMPES-OpenFOAM, nCells = '+str(len(x))+r', ExecutionTime$\approx$'+str(ExeTime_OF) + ' min' + '\nFOSSIL, nCells = '+str(int(len(x_FOSSIL)/4))+r', ExecutionTime$\approx$'+str(ExeTime_FOSSIL)+' min')   

    fig.tight_layout()

    plt.savefig(f'PVI_{tD}.png', dpi=500)
    plt.savefig(f'PVI_{tD}.pdf', dpi=300)
    plt.show()
