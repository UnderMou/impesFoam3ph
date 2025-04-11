from lmfit import Minimizer, Parameters, create_params, report_fit
import os
from jinja2 import Environment, FileSystemLoader
import subprocess
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from aux_funcs import *

def S_norm(S,Sr,Swc,Sgr,Sor):
    return np.divide(S-Sr,1-Swc-Sor-Sgr)

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
                    a_exp = None,
                    b_exp = None,
                    c_exp = None,
                    kra_max = None,
                    krb_max = None, 
                    krc_max = None,
                    swc = None,
                    sgr = None,
                    sor = None,
                   verbose=True):
    """
    Runs an OpenFOAM simulation with the given parameters.

    Parameters:
        base_dir (str): Path to the base OpenFOAM case directory.
        experiment_name (str): Name of the experiment to create.
    """

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

    output = template.render(a_exp = a_exp,
                             b_exp = b_exp,
                             c_exp = c_exp,
                             kra_max = kra_max,
                             krb_max = krb_max,
                             krc_max = krc_max,
                             swc = swc,
                             sgr = sgr,
                             sor = sor
                            )

    # Overwrite the template transportProperties file with the new values
    with open(os.path.join(new_dir, 'constant', 'transportProperties'), 'w') as f:
        f.write(output)

    try:
        print(os.getcwd())
        result = subprocess.run(
            ["bash", "run_solver.sh", new_dir],  
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error running the script:", e.stderr)
    
def residual(params, base_case_directory, experiment_name, idt, CTdata):
    parvals = params.valuesdict()
    a_exp = parvals['a_exp']
    b_exp = parvals['b_exp']
    c_exp = parvals['c_exp']
    kra_max = parvals['kra_max']
    krb_max = parvals['krb_max']
    krc_max = parvals['krc_max']
    swc = parvals['swc']
    sgr = parvals['sgr']
    sor = parvals['sor']

    print('run!')

    run_simulation(base_case_directory, experiment_name, 
                   a_exp = a_exp,
                   b_exp = b_exp,
                   c_exp = c_exp,
                   kra_max = kra_max,
                   krb_max = krb_max, 
                   krc_max = krc_max,
                   swc = swc,
                   sgr = sgr,
                   sor = sor,
                   verbose=True)
    
    nGet = -2
    res = []
    for i in range(len(idt)):
        data_dict = parse_openfoam_case('../simulator/experiments/'+experiment_name, variables=['Sa','Sb'])
        Sa = data_dict['Sa'].iloc[idt[i]][:nGet]
        Sb = data_dict['Sb'].iloc[idt[i]][:nGet]
        Sc = 1.0 - Sa - Sb

        Sa = S_norm(Sa,sgr,swc,sgr,sor)
        Sb = S_norm(Sb,swc,swc,sgr,sor)
        Sc = S_norm(Sc,sor,swc,sgr,sor)

        x = np.linspace(0,1,len(Sa))
        L = x[-1]
        x = x/L
        idx = np.argmin(np.abs(x-0.5))

        SoCT = np.interp(x,CTdata[i]['x_SoCT'],CTdata[i]['So_CT'])
        SgCT = np.interp(x,CTdata[i]['x_SgCT'],CTdata[i]['Sg_CT'])
        SwCT = np.interp(x,CTdata[i]['x_SwCT'],CTdata[i]['Sw_CT'])

        error_Sw = np.abs(SwCT - Sb)/np.max(SwCT)
        error_Sg = np.abs(SgCT - Sa)/np.max(SgCT)
        error_So = np.power(SoCT - Sc, 2)
        error_So = error_So[:idx]

        plt.figure()
        plt.plot(x,SoCT,'b')
        plt.scatter(CTdata[i]['x_SoCT'],CTdata[i]['So_CT'],c='k')
        plt.plot(x,Sc,'gray')
        error = error_So
        plt.scatter(x[:idx],error,c='r')
        plt.savefig('plot.png',dpi=200)
        # plt.show()
        plt.close()

        res.append(np.concatenate((error_Sg, error_So, error_Sw)))
    
    # print(np.array(res).shape)
    # print(np.sum(res))
    # exit()
    print(np.sum(error_So))

    return error_So # res

if __name__ == "__main__":

    experiment_name = 'fitting_lmfit'
    BaseCase_dir = '/home/anderson/OpenFOAM/anderson-9/run/uqsa/impesFoam3ph/simulator/base_cases/Lyu_CT_exp_lmfit/'
    new_dir = "../simulator/experiments/" + experiment_name

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    else:
        print(" -- The directory already exists. Files will be overwritten. --")

    # Time inputs
    desired_times = [800]
    ti = 0
    tf = 2000
    write_interval = 200
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    idt = [np.argmin(np.abs(t-t_find)) for t_find in desired_times]

    # CT Scan data
    tD = [0.36]
    CTscan_results = []
    for i in range(len(tD)):
        Sg_CT_file = 'Tang/Sg_PVI_'+str(tD[i])+'.csv'
        Sw_CT_file = 'Tang/Sw_PVI_'+str(tD[i])+'.csv'
        So_CT_file = 'Tang/So_PVI_'+str(tD[i])+'.csv'
        CTscan_results.append(get_CTscanData(So_CT_file, Sg_CT_file, Sw_CT_file))

    # Parameters ranges
    a_expR      =  [1,5] # 1.62 
    b_expR      =  [1,5] # 3.86 
    c_expR      =  [1,5] # 2.54 
    kra_maxR    =  [0.1,1] # 0.83 
    krb_maxR    =  [0.1,1] # 0.247 
    krc_maxR    =  [0.1,1] # 0.584 
    swcR        =  [0.01,0.2] # 0.197 
    sgrR        =  [0.01,0.2] # 0.013 
    sorR        =  [0.01,0.2] # 0.103 
    range_val   = [a_expR,b_expR,c_expR,kra_maxR,krb_maxR,krc_maxR,swcR,sgrR,sorR]   
    
    # Parameters initial
    a_expI      =  1.917 # 1.62 
    b_expI      =  3.451 # 3.86 
    c_expI      =  2.899 # 2.54 
    kra_maxI    =  0.757 # 0.83 
    krb_maxI    =  0.293 # 0.247 
    krc_maxI    =  0.532 # 0.584 
    swcI        =  0.197 
    sgrI        =  0.013 
    sorI        =  0.103 
    init_val    = [a_expI,b_expI,c_expI,kra_maxI,krb_maxI,krc_maxI,swcI,sgrI,sorI]   

    # Parameters search space
    param_labels = ['a_exp', 'b_exp', 'c_exp', 'kra_max', 'krb_max', 'krc_max', 'swc', 'sgr', 'sor']
    params = Parameters()
    for i in range(len(param_labels)):
        if param_labels[i] in ['swc', 'sgr', 'sor']:
            params.add(param_labels[i], vary=False, value=init_val[i],min=range_val[i][0] , max=range_val[i][1])
        else:
            params.add(param_labels[i], vary=True, value=init_val[i],min=range_val[i][0] , max=range_val[i][1])
    print(params)

    # Parameter fitting 
    foo = Minimizer(residual, params, fcn_kws={'base_case_directory': BaseCase_dir, 'experiment_name': experiment_name, 'idt': idt, 'CTdata': CTscan_results})
    result = foo.minimize(method='differential_evolution')

    # residual(params, base_case_directory= BaseCase_dir, experiment_name= experiment_name)