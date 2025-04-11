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
from uqpylab import sessions, display_util

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
                   mu_a = None,
                   mu_b = None,
                   mu_c = None,
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

    output = template.render(mu_a = mu_a, 
                             mu_b = mu_b,
                             mu_c = mu_c,
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
    i, (mu_a, mu_b, mu_c, fmoil, floil, epoil) = params
    experiment_name = experiment_name + f"/sample_{i:03d}"
    run_simulation(base_case_directory, experiment_name, 
                    mu_a = mu_a, 
                    mu_b = mu_b,
                    mu_c = mu_c,
                    fmoil=fmoil,
                    floil=floil,
                    epoil=epoil,
                    verbose=False)

if __name__ == '__main__':

    Swc = 0.197
    Sgr = 0.013
    Sro = 0.103

    # INPUT DIRECTORIES
    experiment_name = 'SA_CTscan_ViscoFoil'
    BaseCase_dir = '/home/anderson/OpenFOAM/anderson-9/run/uqsa/impesFoam3ph/simulator/base_cases/Lyu_CT_exp_ViscoFoil/'
    new_dir = "../simulator/experiments/" + experiment_name

    ############ UQ[py]Lab - INPUT #############

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
    mu_aR = 1.5e-5
    mu_bR = 0.7e-3
    mu_cR = 9.0e-2
    fmoilR = 0.823
    delta = 0.528
    epoilR = 3.827

    dpr = 0.2

    dists = {
    'num_vars': 6,
    'names': ['mu_g', 'mu_w', 'mu_o', 'fmoil', 'delta', 'epoil'],
    'bounds': [[mu_aR*(1-dpr), mu_aR*(1+dpr)],
               [mu_bR*(1-dpr), mu_bR*(1+dpr)],
               [mu_cR*(1-dpr), mu_cR*(1+dpr)],
               [fmoilR*(1-dpr), fmoilR*(1+dpr)],
               [delta*(1-dpr), delta*(1+dpr)],
               [epoilR*(1-dpr), epoilR*(1+dpr)]
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
            }
        ]
    }

    myInput = uq.createInput(InputOpts)

    uq.print(myInput)
    # uq.display(myInput)

    nSamples = 1000
    X_ED = uq.getSample(myInput,nSamples)
    print(type(X_ED), "\n", X_ED)

    mu_aR_samples = np.array([sample[0] for sample in X_ED])
    mu_bR_samples = np.array([sample[1] for sample in X_ED])
    mu_cR_samples = np.array([sample[2] for sample in X_ED]) 
    fmoil_samples = np.array([sample[3] for sample in X_ED])
    delta_samples = np.array([sample[4] for sample in X_ED])
    epoil_samples = np.array([sample[5] for sample in X_ED]) 
    
    # DEFINE floil as fmoil > floil
    floil_samples = fmoil_samples * (1.0 - delta_samples)
    check = floil_samples > fmoil_samples
    print(check)
    print("is any floil > fmoil: ", any(check))
    X_ED[:,4] = floil_samples

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    else:
        print(" -- The directory already exists. Files will be overwritten. --")

    np.savetxt(new_dir + '/X_ED.csv',X_ED,delimiter=',')

    data = pd.DataFrame({
        'mu_aR': mu_aR_samples,
        'mu_bR': mu_bR_samples,
        'mu_cR': mu_cR_samples,
        'fmoil': fmoil_samples,
        'floil': floil_samples,
        'epoil': epoil_samples
    })
    
    sns.pairplot(data)
    plt.savefig(new_dir + '/pairplot.pdf', dpi=200)
    plt.show()
       
    process_func = partial(
    process_simulation, 
    base_case_directory=BaseCase_dir,
    experiment_name=experiment_name
    )

    params = list(enumerate(X_ED))

    # Continuing simulation
    run_from = 136
    params = params[run_from:]
    print(params)

    nthreads = 6

    with Pool(nthreads) as pool:
        for _ in tqdm(
            pool.imap_unordered(process_func, params),
            total=len(params), 
            desc='Running simulations',
            mininterval=1.0     # Updates at most once per second
        ):
            pass