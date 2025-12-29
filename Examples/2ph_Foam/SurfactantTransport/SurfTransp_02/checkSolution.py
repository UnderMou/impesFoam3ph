import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import scienceplots
import os

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
    # data.index = data.index.astype(int)

    return data

if __name__ == '__main__':

    # Parameters
    Swc = 0.437
    Sgr = 0.293
    phi = 0.21
    L = 1.0
    Cmax = 2


    plt.style.use('science')
    plt.rcParams.update({'font.size': 18})

    data_sw = pd.read_csv('Sw_analytic_01.csv',header=None)
    data_cs = pd.read_csv('Cs_analytic_01.csv',header=None)

    time_dirs = np.linspace(0,0.02835,int(0.02835/0.001)+1).tolist()
    time_dirs = [round(t, 3) for t in time_dirs]
    time_dirs[0] = 0

    print(time_dirs[1])

    data_dict = parse_openfoam_case('', variables=['Sb','Cs'],time_dirs=time_dirs)

    Sb = data_dict['Sb'].iloc[-1]
    Sb_norm = (Sb - Swc) / (1 - Swc - Sgr)
    Cs = data_dict['Cs'].iloc[-1]
    Cs_norm = Cs/Cmax
    x = np.linspace(0,1,len(Sb))

    fig, ax = plt.subplots(1,2,figsize=(10,5))

    ax[0].plot(data_sw.iloc[:,0],data_sw.iloc[:,1], lw=2, c='b', label='Fritis et al. (2023)')
    ax[0].scatter(x, Sb_norm, c='k', label='ImpesFOAM')
    ax[0].set_xlabel(r'$x_D$ [-]')
    ax[0].set_ylabel(r'$S$ [-]')
    ax[0].set_ylim([0,1])
    ax[0].grid(True)
    ax[0].legend(fontsize=14, loc='lower right')

    ax[1].plot(data_cs.iloc[:,0],data_cs.iloc[:,1], lw=2, c='b')
    ax[1].scatter(x, Cs_norm, c='k')
    ax[1].set_ylim([0,1])
    ax[1].set_xlabel(r'$x_D$ [-]')
    ax[1].set_ylabel(r'$C$ [-]')
    ax[1].set_ylim([0,1])
    ax[1].grid(True)

    plt.tight_layout()
    plt.savefig('surfactant_02.pdf', dpi=300)
    plt.show()