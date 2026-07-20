import matplotlib
matplotlib.use('Agg')   # sem interface gráfica

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

def parse_openfoam_case(case_dir, variables=['p','K', 'Sa', 'Sb', 'Sc', 'U', 'Ua', 'Ub', 'Uc', 'Fa', 'Fb'], time_dirs=None):
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

    ti = 0
    tf = 1e8
    write_interval = 5e6
    t = np.linspace(ti, tf, int(tf/write_interval)+1)

    nx, ny = 110, 30
    lx, ly = 360, 100

    insts = [2e7, 5e7]

    dx = lx/nx
    dy = ly/ny
    dz = 10.0
    dv = dx*dy*dz
    qt = 4.883307965499746e-08
    Qt = ny * qt * dv

    phi = 0.07
    PV = lx * ly * 10.0 * phi

    fig, axes = plt.subplots(2, 2, figsize=(12, 6), constrained_layout=True)

    dirs = ['fg70_Close2corner', 'fg70_corner']

    cond = ['A', 'B']

    # Plot the oil phase saturation maps

    for i, d in enumerate(dirs):

        data_dict = parse_openfoam_case(d, variables=['p','K','Sc'])
        data_dict.index = (data_dict.index.astype(float) / write_interval).astype(int)

        sats_OF = data_dict['Sc']

        insts_id = np.array([np.argmin(np.abs(t-inst)) for inst in insts], dtype=int)

        pvis = Qt/PV * np.array(insts)

        sat_OF = [np.reshape(sats_OF[j], (ny, nx)) for j in insts_id]

        # loop over columns (time)
        for k in range(2):
            ax = axes[i, k]

            im = ax.imshow(
                sat_OF[k],
                cmap='coolwarm',
                origin='lower',
                vmin=0.12, vmax=0.79,
                extent=[0, lx, 0, ly],
                aspect='auto'
            )

            # --- Labels ---
            if i == 1:  # bottom row
                ax.set_xlabel(r'$x$ [m]')
            else:
                ax.set_xticklabels([])

            if k == 0:  # left column
                ax.set_ylabel(r'$y$ [m]')
            else:
                ax.set_yticklabels([])

            # --- Titles (top row only) ---
            if i == 0:
                ax.set_title(f'{pvis[k]:.2f} PVI',fontsize=16)

        axes[i, 0].text(
            -0.15, 0.5, cond[i],
            transform=axes[i, 0].transAxes,
            rotation=90,
            va='center',
            ha='center',
            fontsize=14
        )

    # Better colorbar
    cbar = fig.colorbar(
        im,
        ax=axes,
        orientation='horizontal',
        fraction=0.03,
        pad=0.03
    )

    cbar.set_label(r'$S_o$', fontsize=14)

    plt.savefig('soi_heterog_spe.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    # Plot the Gd map

    mu_o = 1.1e-3   # [Pa.s]
    rho_w = 1028.0  # [kg/m3]
    rho_o = 750.0   # [kg/m3]
    delta_rho = rho_w - rho_o

    fig, axes = plt.subplots(2, 2, figsize=(12, 6), constrained_layout=True)

    for i, d in enumerate(dirs):

        data_dict = parse_openfoam_case(d, variables=['p','K','Sc','U'])
        data_dict.index = (data_dict.index.astype(float) / write_interval).astype(int)

        K_OF = data_dict['K']
        U_OF = data_dict['U']

        insts_id = np.array([np.argmin(np.abs(t-inst)) for inst in insts], dtype=int)

        pvis = Qt/PV * np.array(insts)

        K_OF = np.array([np.reshape(K_OF[j], (ny, nx)) for j in insts_id])
        U_OF = np.array([np.reshape(U_OF[j], (ny, nx, 3)) for j in insts_id])
        Umag_OF = np.linalg.norm(U_OF, axis=-1)

        # print(K_OF.shape)
        # exit()
    

        Gd_OF = K_OF * delta_rho * 9.81 / (mu_o * Umag_OF)
   
        # loop over columns (time)
        for k in range(2):
            ax = axes[i, k]

            im = ax.imshow(
                Gd_OF[k],
                cmap='coolwarm',
                origin='lower',
                vmin=0.01, vmax=121.86,
                extent=[0, lx, 0, ly],
                aspect='auto'
            )

            # --- Labels ---
            if i == 2:  # bottom row
                ax.set_xlabel(r'$x$ [m]')
            else:
                ax.set_xticklabels([])

            if k == 0:  # left column
                ax.set_ylabel(r'$y$ [m]')
            else:
                ax.set_yticklabels([])

            # --- Titles (top row only) ---
            if i == 0:
                ax.set_title(f'{pvis[k]:.2f} PVI',fontsize=16)

        axes[i, 0].text(
            -0.15, 0.5, cond[i],
            transform=axes[i, 0].transAxes,
            rotation=90,
            va='center',
            ha='center',
            fontsize=14
        )

    # Better colorbar
    cbar = fig.colorbar(
        im,
        ax=axes,
        orientation='horizontal',
        fraction=0.03,
        pad=0.03
    )

    cbar.set_label(r'$G_d$', fontsize=14)

    plt.savefig('Gd_heterog_spe.pdf', dpi=300, bbox_inches='tight')
    plt.show()