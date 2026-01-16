import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
from matplotlib.lines import Line2D

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

def S_norm(S,Sr,Swc,Sgr,Sor):
    return np.divide(S-Sr,1-Swc-Sor-Sgr)

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams['font.size'] = 16

    # Reading Analytical solution from Lyu
    data_Lyu_Sg = pd.read_csv('Sg_analytical_LyuC3.2.csv')
    data_Lyu_So = pd.read_csv('So_analytical_LyuC3.2.csv')

    # Define t vector
    ti = 0
    tf = 60000
    write_interval = 50
    t = np.linspace(ti,tf,int(tf/write_interval)+1)
    print(t)
    t[-1] = int(1)
    time_dirs = [f'{x:.3f}'.rstrip('0').rstrip('.') if x % 1 else f'{int(x)}' for x in t]
    print(time_dirs)

    path = ''
    data_dict = parse_openfoam_case(path, variables=['Sb','Sa'], time_dirs=time_dirs)

    t_find = 53500
    id_t = np.argmin(np.abs(t-t_find))
    print(t[id_t])
    Sw_OF = data_dict['Sb'].iloc[id_t]
    Sg_OF = data_dict['Sa'].iloc[id_t]
    So_OF = 1.0 - Sw_OF - Sg_OF
    

    nAvoid = -1
    Sw_OF = Sw_OF[:nAvoid]
    Sg_OF = Sg_OF[:nAvoid]
    So_OF = So_OF[:nAvoid]
    xD = np.linspace(0,1,len(Sw_OF))

    Swc = 0.1
    Sgr = 0
    Sor = 0.1
    Sg_OF = S_norm(Sg_OF,Sgr,Swc,Sgr,Sor)
    So_OF = S_norm(So_OF,Sor,Swc,Sgr,Sor)
    
    # Plotting the results
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5))

    # ax1.plot(xD_Fo, Sg_Fo, linestyle='--', color='b', label='FOSSIL')
    ax1.plot(xD, Sg_OF, lw=2,linestyle='--', color='k', label='ImpesFOAM')
    ax1.plot(data_Lyu_Sg.iloc[:,0], data_Lyu_Sg.iloc[:,1], color='g', label='Analytical - Tang et al., 2022')
    # ax1.set_xlabel(r'$x_D$')
    ax1.set_ylabel(r'$S_g$')
    # ax1.set_title(r'$S_g$ profile at $T_D=1.0$')
    ax1.grid()
    ax1.legend(fontsize=18)
    # ax1.set_ylim([-0.05,0.85])

    # ax2.plot(xD_Fo, So_Fo, linestyle='--', color='b', label='FOSSIL')
    ax2.plot(xD, So_OF, lw=2,linestyle='--', color='k', label='ImpesFOAM')
    ax2.plot(data_Lyu_So.iloc[:,0], data_Lyu_So.iloc[:,1], color='g')
    # ax2.set_xlabel(r'$x_D$')
    ax2.set_ylabel(r'$S_o$')
    # ax2.set_title(r'$S_o$ profile at $T_D=1.0$')
    ax2.grid()
    # ax2.set_ylim([-0.05,0.85])

    plt.tight_layout()
    plt.savefig('LyuC3.2.pdf', dpi=300)
    plt.show()
