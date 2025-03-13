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

if __name__ == '__main__':
    plt.style.use('science')
    plt.rcParams.update({"font.size":15})  

    ti = 0
    tf = 139000
    write_interval = 1000
    t = np.linspace(ti,tf,int(tf/write_interval + 1))

    dir = '.'
    data_dict = parse_openfoam_case(dir, variables=['Fa', 'p'])
    print(data_dict.head())
    get = -2
    press_drop = data_dict['p'].apply(lambda x: x[0] - x[get])
    press_drop = press_drop.to_numpy()
    print(press_drop)
    time_steps = len(t)-1
    t = t[:-1]

    plt.figure(figsize=(10,7))
    plt.plot(t,press_drop)
    plt.ylabel(r'$\Delta p$ - pressure drop [Pa]')
    plt.xlabel(r'time [s]')
    plt.grid()
    plt.show()
    