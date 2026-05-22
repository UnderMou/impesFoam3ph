import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from jinja2 import Environment, FileSystemLoader

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

def run_simulation(base_dir, experiment_name, mu_a=1.0e-5, mu_b=1.0e-3, mu_c=5.0e-2, verbose=True):
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

    output = template.render(mu_a=mu_a, mu_b=mu_b, mu_c=mu_c)

    # Overwrite the template transportProperties file with the new values
    with open(os.path.join(new_dir, 'constant', 'transportProperties'), 'w') as f:
        f.write(output)

    # try:
    #     result = subprocess.run(
    #         ["bash", "run_solver.sh", new_dir],  
    #         check=True,
    #         capture_output=True,
    #         text=True
    #     )
    #     print(result.stdout)
    # except subprocess.CalledProcessError as e:
    #     print("Error running the script:", e.stderr)

if __name__ == '__main__':

    # sup_data_dir = os.getcwd() + '/sup_data.csv'
    # sup_data = pd.read_csv(sup_data_dir)
    # # print(sup_data.head())
    # x = sup_data['Points:2'][:-1]
    # print(x.shape)

    # # case_dir = f'/home/anderson/OpenFOAM/anderson-9/run/impesFoam3ph_v1/Lyu_CT_exp'
    # case_dir = os.getcwd()
    # data_dict = parse_openfoam_case(case_dir, variables=['Sa','Sb'])

    # # print(data_dict.shape)
    # # print(list(data_dict.columns.values))

    # ti = 0
    # tf = 10000
    # write_interval = 25
    # t = np.linspace(ti,tf,int(tf/write_interval + 1))
    # print(t)
    # t_find = 1350 # tD = 0.36
    # idx = np.argmin(np.abs(t-t_find))
    # print(idx)

    # Sa = data_dict['Sa'].iloc[idx]
    # Sb = data_dict['Sb'].iloc[idx]
    # Sc = 1.0 - Sa - Sb
    # print(Sa.shape)

    # plt.plot(x,Sc)
    # plt.show()

    case_dir = os.getcwd()
    experiment_name = 'teste'
    run_simulation(case_dir, experiment_name, mu_a=1.0e-5, mu_b=1.0e-3, mu_c=5.0e-2, verbose=True)