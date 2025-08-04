import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scienceplots

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

if __name__ == "__main__":
    plt.style.use('science')
    plt.rcParams.update({'font.size': 16})

    # t_find = 0.5 # 0.25
    t_find = [0.2, 0.4]

    # Domain:
    x = np.linspace(0,1,500)

    # Time controls:
    ti = 0
    tf = 1
    write_interval = 0.01
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    
    
    
    time_dirs = np.linspace(0,1,int(1/0.01)+1).tolist()
    time_dirs = [round(t, 2) for t in time_dirs]
    time_dirs[0] = 0
    time_dirs[-1] = 1
    # print(time_dirs)
    # print(time_dirs[idt])

    labels = [0.001, 0.005, 0.01, 0.02]
    exps = ['eps_0.001',  'eps_0.005',  'eps_0.01',  'eps_0.02']
    style = ['-', '-', '--', ':']
    width = [1,0.5,1,1]

    fig, ax = plt.subplots(1,3,figsize=(15,5))
    for j in range(len(t_find)):
        for i in range(len(exps)):

            data_dict = parse_openfoam_case(exps[i], variables=['Sb','Sa'],time_dirs=time_dirs)

            idt = np.argmin(np.abs(t-t_find[j]))

            Sb = np.array(data_dict['Sb'].iloc[idt])
            Sa = np.array(data_dict['Sa'].iloc[idt])
            Sc = 1 - Sb - Sa

        
            if j == 0: ax[0].plot(x,Sb,c='k',linewidth=width[i],linestyle=style[i],label=f'$\epsilon = {labels[i]}$')
            else: ax[0].plot(x,Sb,c='k',linewidth=width[i],linestyle=style[i])
            ax[1].plot(x,Sa,c='k',linewidth=width[i],linestyle=style[i])
            ax[2].plot(x,Sc,c='k',linewidth=width[i],linestyle=style[i])

    ax[0].set_xlim(0, 1)       
    ax[0].set_ylim(0, 1)        
    ax[0].set_xlabel(r'$x$ [m]')
    ax[0].set_ylabel(r'$S_w$ [-]')
    ax[0].grid(True)

    ax[1].set_xlim(0, 1)       
    ax[1].set_ylim(0, 1)       
    ax[1].set_xlabel(r'$x$ [m]')
    ax[1].set_ylabel(r'$S_g$ [-]')
    ax[1].grid(True)

    ax[2].set_xlim(0, 1)       
    ax[2].set_ylim(0, 1)     
    ax[2].set_xlabel(r'$x$ [m]')
    ax[2].set_ylabel(r'$S_o$ [-]')
    ax[2].grid(True)

    ax[0].legend(loc='upper right')

    plt.tight_layout()
    plt.savefig('capillaryPressRP1.pdf', dpi=300)
    plt.show()