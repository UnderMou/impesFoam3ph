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

    plt.style.use('science')
    plt.rcParams.update({'font.size': 18})

    data_sw = pd.read_csv('Sw_Kc001.csv')
    data_nd = pd.read_csv('nD_Kc001.csv')

    # Shock front velocity
    u = 2.93e-5
    phi = 0.25

    fw_plus = 1.0
    fw_minus = 0.268
    Sw_plus = 0.72
    Sw_minus = 0.372
    nD_minus = 0.664
    nD_plus = 1.0
    vs = (u/phi)*(fw_plus-fw_minus)/(Sw_plus-Sw_minus)                    

    # Time and domain
    dt = 0.25
    time_dirs = np.arange(dt,70+dt,dt)
    time_dirs = [int(x) if x.is_integer() else float(x) for x in time_dirs]
    print(time_dirs)
    data_dict = parse_openfoam_case('', variables=['Sb','nD'],time_dirs=time_dirs)

    t_find = 150
    idt = np.argmin(np.abs(np.array(time_dirs) - t_find))

    x = np.linspace(0,0.051,len(data_dict['Sb'][0]))
    x = x[:-2]
    t_calc = 58
    eta = (x - vs*t_calc)*1000 # m to mm

    fig, ax = plt.subplots(2,1,figsize=(10,5))

    Sw = data_dict['Sb'].iloc[idt]
    Sw = Sw[:-2]
    nD = data_dict['nD'].iloc[idt]
    nD = nD[:-2]
    

    ax[0].plot(eta,Sw, lw=2, c='k',label='ImpesFOAM')   
    ax[0].plot(np.concatenate([[eta[0]],(data_sw.iloc[:,0]),[eta[-1]]]),np.concatenate([[Sw_minus],data_sw.iloc[:,1],[Sw_plus]]),linestyle='--', lw=2, c='r',label='Ashoori et al. (2011)')      
    ax[0].set_xlabel(r'$\eta$ [mm]')
    ax[0].set_ylabel(r'$S_w$ [-]')
    # ax[0].set_xlim([-3,3])
    ax[0].set_ylim([0,1])
    ax[0].grid(True)
    ax[0].legend(fontsize=14)

    ax[1].plot(eta,nD, lw=2, c='k')     
    ax[1].plot(np.concatenate([[eta[0]],(data_nd.iloc[:,0]),[eta[-1]]]),np.concatenate([[nD_minus],data_nd.iloc[:,1],[nD_plus]]),linestyle='--', lw=2, c='r')       
    ax[1].set_xlabel(r'$\eta$ [mm]')
    ax[1].set_ylabel(r'$n_D$ [-]')
    # ax[1].set_xlim([-3,3])
    ax[1].set_ylim([0.6,1.01])
    ax[1].grid(True)

    plt.legend()
    plt.tight_layout()
    plt.savefig('Ashoori_Kc0_01.pdf', dpi=300)
    plt.show()