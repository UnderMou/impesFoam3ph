import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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
        typeOf = content[start_index].split()[1]    

        # Check if the field was written as uniform - special handling
        if typeOf == 'uniform':
            num_elements = 100 # TODO: Generalizar para outras malhas 
            value = float(content[start_index].split()[2][0])
            values = value * np.ones(num_elements)
            return values
        elif typeOf == 'nonuniform':
            num_elements = int(content[start_index + 1]) 
        else:
            raise ValueError("Check further")
        
        
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
        exit()
        return None

def parse_openfoam_case(case_dir, variables=['p', 'Sa', 'Sb', 'U', 'Ua', 'Ub', 'Uc', 'Fa', 'Fb','Fshear','Foil'], time_dirs=None):
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

if __name__ == "__main__":

    t_find = 250

    # Parameters
    rho_w = 1e3
    rho_g = 1.0
    g = 9.81

    # Domain:
    x = np.linspace(0,1,100)

    # Time controls:
    ti = 0
    tf = 2e6 
    write_interval = 50
    t = np.linspace(ti,tf,int(tf/write_interval + 1))
    idt = np.argmin(np.abs(t-t_find))

    data_dict = parse_openfoam_case('.', variables=['Sb','dpcdS'])
    Sb = np.array([data_dict['Sb'].iloc[i] for i in range(idt)])
    Sa = np.array([1-sb for sb in Sb])
    dpcdS = np.array([data_dict['dpcdS'].iloc[i] for i in range(idt)])
    dSdx = np.array([((rho_g-rho_w)*g)/dpcdS_ for dpcdS_ in dpcdS])

    analytic_data = pd.read_csv('dSdx_analytic.csv', header=None)
    an_S = analytic_data.iloc[:,0]
    an_dSdx = analytic_data.iloc[:,1]

    analytic_data = pd.read_csv('Sw_analytic.csv', header=None)
    an_x_prof = analytic_data.iloc[:,0]
    an_S_prof = analytic_data.iloc[:,1]
    
    S = np.linspace(np.min(Sb), np.max(Sb), len(x))

    # Set up the figure and axis
    fig, ax = plt.subplots(1,2, figsize=(7,5))
    line_Sw = ax[0].scatter([], [], label='$S_w$')
    line_Sg = ax[0].scatter([], [], label='$S_g$')
    line_dSdx = ax[1].scatter([], [], label='$dS_w/dx$')
    line_andSdx, = ax[1].plot([], [], color='k', lw=2, label='analytical $dS_w/dx$')
    line_anSw, = ax[0].plot([], [], color='k', lw=2, label='steady-state $S_w(x)$\nCarrillo et al. (2020)')
    ax[0].set_xlim(-.1, 1.1)        # ax.set_xlim(x[0], x[-1])
    ax[0].set_ylim(-.1, 1.1)        # ax.set_ylim(Sb.min(), Sb.max())
    ax[0].set_xlabel(r'$x$ [m]')
    ax[0].set_ylabel(r'$S$ [-]')
    ax[0].grid(True)
    ax[0].legend(loc='center right')

    ax[1].set_xlim(-.1, 1.1)        # ax.set_xlim(x[0], x[-1])
    ax[1].set_ylim(-.1, 25)         # ax.set_ylim(Sb.min(), Sb.max())
    ax[1].set_xlabel(r'$S_w$ [-]')
    ax[1].set_ylabel(r'$dS_w/dx$ [m$^{-1}$]')
    ax[1].grid(True)
    ax[1].legend(loc='upper left')

    plt.tight_layout()

    # Init function: clears the line
    def init():
        line_Sw.set_offsets(np.empty((0, 2)))
        line_Sg.set_offsets(np.empty((0, 2)))
        line_dSdx.set_offsets(np.empty((0, 2)))
        line_andSdx.set_data([], [])
        line_anSw.set_data([], [])
        return line_Sw, line_Sg, line_dSdx, line_andSdx, line_anSw,

    # Animation function: updates line
    def animate(i):
        line_Sw.set_offsets(np.column_stack((x, Sb[i])))
        line_Sg.set_offsets(np.column_stack((x, Sa[i])))
        line_dSdx.set_offsets(np.column_stack((Sb[i], dSdx[i])))
        line_andSdx.set_data(an_S, an_dSdx)
        line_anSw.set_data(an_x_prof, an_S_prof)
        ax[0].set_title(f"Time step {i}")
        ax[1].set_title(f"Time step {i}")
        
        return line_Sw,line_Sg,line_dSdx,line_andSdx,line_anSw,

    # Create animation
    ani = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=idt, interval=100, blit=True)

    plt.show()
    plt.close()

    # Steady-state plot
    plt.rcParams.update({'font.size': 18})
    fig2, ax2 = plt.subplots(1,2, figsize=(8,5))

    ax2[0].scatter(x, Sb[-1], label='$S_w$')
    ax2[0].scatter(x, Sa[-1], label='$S_g$')
    ax2[0].plot(an_x_prof, an_S_prof, color='k', lw=2, label='steady-state $S_w(x)$\nCarrillo et al. (2020)')

    ax2[1].scatter(an_S, an_dSdx, label='$dS_w/dx$')
    ax2[1].plot(an_S, an_dSdx, color='k', lw=2, label='analytical $dS_w/dx$')

    ax2[0].set_xlim(-.1, 1.1)        # ax.set_xlim(x[0], x[-1])
    ax2[0].set_ylim(-.1, 1.1)        # ax.set_ylim(Sb.min(), Sb.max())
    ax2[0].set_xlabel(r'$x$ [m]')
    ax2[0].set_ylabel(r'$S$ [-]')
    ax2[0].grid(True)
    ax2[0].legend(loc='center right', fontsize=10)

    ax2[1].set_xlim(-.1, 1.1)        # ax.set_xlim(x[0], x[-1])
    ax2[1].set_ylim(-.1, 25)         # ax.set_ylim(Sb.min(), Sb.max())
    ax2[1].set_xlabel(r'$S_w$ [-]')
    ax2[1].set_ylabel(r'$dS_w/dx$ [m$^{-1}$]')
    ax2[1].grid(True)
    ax2[1].legend(loc='upper left', fontsize=10)
    plt.tight_layout()
    plt.show()