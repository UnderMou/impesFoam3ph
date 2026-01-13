import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import scienceplots
import matplotlib.gridspec as gridspec

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
    plt.rcParams.update({'font.size': 14})

    data_dict = parse_openfoam_case('.', variables=['p','Sb'])
    
    ti = 0
    tf = 5e8
    write_interval = 1e5
    t = np.linspace(ti,tf,int(tf/write_interval)+1)
    

    cols = nx = 51
    lines = ny = 51
    lx = 142.245
    ly = 142.245
    lz = 10

    for i in data_dict.index:
        data_dict.rename(index={i: int(i/write_interval)}, inplace=True)


    sats_OF = data_dict['Sb']

    insts = [1e7, 4e7, 7e7, 1.0e8]
    insts_id = np.array([np.argmin(np.abs(t-inst)) for inst in insts], dtype=int)


    sat_OF = [np.reshape(sats_OF[inst_id], (ny, nx)) for inst_id in insts_id]

    # fig, axs = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(7,5))
    # # fig.subplots_adjust(hspace=0.02, wspace=0.02)

    # im1 = plt.imshow(sat_OF, cmap='coolwarm', origin='lower', vmin=0.3, vmax=0.93, extent=[0, lx, 0, ly])
    # plt.xlabel(r'$x$ [m]')
    # plt.ylabel(r'$y$ [m]')

    # # Configura o espaço e a posição da colorbar ao lado direito da figura
    # cbar_ax = fig.add_axes([0.87, 0.35, 0.02, 0.418])
    # fig.colorbar(im1, cax=cbar_ax, label='$S_w$')
    # cbar_ax.yaxis.label.set_size(16)
    # cbar_ax.tick_params(labelsize=16)
    
    # plt.tight_layout()
    # plt.show()




    # Cria figura e especifica layout com GridSpec
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 0.05], hspace=0.4, wspace=0.2)

    # Subplots 2x2 nas duas primeiras linhas
    axes = []
    aux = 0
    for i in range(2):
        for j in range(2):
            ax = fig.add_subplot(gs[i, j])
            im = ax.imshow(sat_OF[aux], cmap='coolwarm', origin='lower',
                        vmin=0.3, vmax=0.93, extent=[0, lx, 0, ly])
            ax.set_xlabel(r'$x$ [m]', labelpad=6)
            ax.set_ylabel(r'$y$ [m]', labelpad=6)
            axes.append(ax)
            aux += 1

    # Colorbar ocupando a linha 3 (horizontal)
    cbar_ax = fig.add_subplot(gs[2, :])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('$S_w$', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig('5spots_4timesPlot.pdf', dpi=300)
    plt.show()