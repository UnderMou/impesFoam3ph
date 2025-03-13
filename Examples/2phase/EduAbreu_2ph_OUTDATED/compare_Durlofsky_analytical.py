import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
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

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams['font.size'] = 14

    # Reading Durlofsky data
    Durlof_Grav_t025 = pd.read_csv('Durlofsky_Grav_t025.csv')
    Durlof_Grav_t050 = pd.read_csv('Durlofsky_Grav_t050.csv')
    Durlof_noGrav_t025 = pd.read_csv('Durlofsky_noGrav_t025.csv')
    Durlof_noGrav_t050 = pd.read_csv('Durlofsky_noGrav_t050.csv')
    
    # Analytical solution
    x_analytical = np.loadtxt("analytical_x.csv", delimiter=",")
    Sw_analitycal_noGrav_t025 = np.loadtxt("Sw_analitycal_noGrav_t025.csv", delimiter=",")
    Sw_analitycal_noGrav_t050 = np.loadtxt("Sw_analitycal_noGrav_t050.csv", delimiter=",")
    Sw_analitycal_Grav_t025 = np.loadtxt("Sw_analitycal_Grav_t025.csv", delimiter=",")
    Sw_analitycal_Grav_t050 = np.loadtxt("Sw_analitycal_Grav_t050.csv", delimiter=",")

    # Define t vector
    ti = 250/1e5
    tf = 12500/1e5
    write_interval = 250/1e5
    t = np.linspace(ti,tf,int(tf/write_interval))
    # print(t)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    #############################################
    #    RESULTS WITH GRAVITY
    #############################################

    # Reading OpenFOAM files
    # tD = 0.25
    path_grav = os.getcwd()
    time_dirs = [0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06,0.0625,0.065,0.0675,0.07,0.0725,0.075,0.0775,0.08,0.0825,0.085,0.0875,0.09,0.0925,0.095,0.0975,0.1,0.1025,0.105,0.1075,0.11,0.1125,0.115,0.1175,0.12,0.1225,0.125]
    data_dict_grav = parse_openfoam_case(path_grav, variables=['Sb'],time_dirs=time_dirs)
    # print(data_dict_grav)
    t_find = 6250/1e5
    id_t = np.argmin(np.abs(t-t_find))
    Sw_t025_grav = data_dict_grav['Sb'].iloc[id_t]
    # tD = 0.5
    t_find = 12500/1e5
    id_t = np.argmin(np.abs(t-t_find))
    Sw_t050_grav = data_dict_grav['Sb'].iloc[id_t]

   

    # x domain 
    x = np.linspace(0,1,len(Sw_t025_grav))
    
    # Plot profiles
    colors = ['b', 'r']
    axes[0].plot(x, Sw_t025_grav, c=colors[0], label=r'impesFoam2ph_G - $t_D = 0.25$')
    axes[0].scatter(Durlof_Grav_t025.iloc[:,0], Durlof_Grav_t025.iloc[:,1], c=colors[0], label=r'DURLOFSKY (1993) - $t_D = 0.25$')
    axes[0].plot(x_analytical, Sw_analitycal_Grav_t025, c=colors[0], linestyle='--', label=r'analytical - $t_D = 0.25$')

    axes[0].plot(x, Sw_t050_grav, c=colors[1], label=r'impesFoam2ph_G - $t_D = 0.5$')
    axes[0].scatter(Durlof_Grav_t050.iloc[:,0], Durlof_Grav_t050.iloc[:,1], c=colors[1], label=r'DURLOFSKY (1993) - $t_D = 0.5$')
    axes[0].plot(x_analytical, Sw_analitycal_Grav_t050, c=colors[1], linestyle='--', label=r'analytical - $t_D = 0.5$')

    axes[0].set_ylim([0,1.4])
    axes[0].set_xlabel('x [m]')
    axes[0].set_ylabel(r'$S_w$ [-]')
    axes[0].grid()
    axes[0].set_title('gravity effects')
    # axes[0].legend(loc='upper right')

    #############################################
    #    RESULTS WITHOUT GRAVITY
    #############################################

    # Reading OpenFOAM files
    # tD = 0.25
    path_nograv = os.getcwd() + '/../EduAbreu_2ph_noG'
    data_dict_nograv = parse_openfoam_case(path_nograv, variables=['Sb'], time_dirs=time_dirs)
    t_find = 6250/1e5
    id_t = np.argmin(np.abs(t-t_find))
    Sw_t025_nograv = data_dict_nograv['Sb'].iloc[id_t]
    # tD = 0.5
    t_find = 12500/1e5
    id_t = np.argmin(np.abs(t-t_find))
    Sw_t050_nograv = data_dict_nograv['Sb'].iloc[id_t]

    # x domain 
    x = np.linspace(0,1,len(Sw_t025_nograv))
    
    # Plot profiles
    colors = ['b', 'r']
    axes[1].plot(x, Sw_t025_nograv, c=colors[0], label=r'impesFoam2ph_G - $t_D = 0.25$')
    axes[1].scatter(Durlof_noGrav_t025.iloc[:,0], Durlof_noGrav_t025.iloc[:,1], c=colors[0], label=r'DURLOFSKY (1993) - $t_D = 0.25$')
    axes[1].plot(x_analytical, Sw_analitycal_noGrav_t025, c=colors[0], linestyle='--', label=r'analytical - $t_D = 0.25$')

    axes[1].plot(x, Sw_t050_nograv, c=colors[1], label=r'impesFoam2ph_G - $t_D = 0.5$')
    axes[1].scatter(Durlof_noGrav_t050.iloc[:,0], Durlof_noGrav_t050.iloc[:,1], c=colors[1], label=r'DURLOFSKY (1993) - $t_D = 0.5$')
    axes[1].plot(x_analytical, Sw_analitycal_noGrav_t050, c=colors[1], linestyle='--', label=r'analytical - $t_D = 0.5$')

    axes[1].set_ylim([0,1.4])
    axes[1].set_xlabel('x [m]')
    axes[1].set_ylabel(r'$S_w$ [-]')
    axes[1].grid()
    axes[1].set_title('no gravity effects')
    axes[1].legend()

    plt.savefig('gravityComparison.pdf', dpi=300)
    plt.show()