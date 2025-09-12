import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import scienceplots
import matplotlib.gridspec as gridspec
import ternary
import matplotlib.tri as tri

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

    path = 'foam_noG'
    # path = 'noFoam_noG'
    data_dict = parse_openfoam_case(path+'/', variables=['p','Sb','Sa'])
    
    ti = 0
    tf = 1.5e8
    write_interval = 1e5
    t = np.linspace(ti,tf,int(tf/write_interval)+1)
    

    cols = nx = 110
    lines = ny = 30
    lx = 360
    ly = 100
    lz = 10

    for i in data_dict.index:
        data_dict.rename(index={i: int(i/write_interval)}, inplace=True)


    Sws_OF_data = data_dict['Sb']
    Sgs_OF_data = data_dict['Sa']

    insts = [1e7, 4e7, 7e7, 1.0e8]
    insts_id = np.array([np.argmin(np.abs(t-inst)) for inst in insts], dtype=int)

    Sws_OF = [np.reshape(Sws_OF_data[inst_id], (ny, nx)) for inst_id in insts_id]
    Sgs_OF = [np.reshape(Sgs_OF_data[inst_id], (ny, nx)) for inst_id in insts_id]
    Sos_OF = [1.0 - Sws_OF[inst_id] - Sgs_OF[inst_id] for inst_id in np.arange(0,len(insts),1)]
    # print(np.array(Sws_OF).shape, np.array(Sgs_OF).shape, np.array(Sos_OF).shape)


    # Define min/max for each channel (could be data min/max or user-defined)
    # # OPTION 1
    # So_min, So_max = 0., 0.44
    # Sg_min, Sg_max = 0., 0.61
    # Sw_min, Sw_max = 0., 0.72
    # OPTION 2
    So_min, So_max = 0., 0.28
    Sg_min, Sg_max = 0., 0.18
    Sw_min, Sw_max = 0., 0.72
    # # OPTION 3
    # So_min, So_max = 0.1, 0.28
    # Sg_min, Sg_max = 0., 0.18
    # Sw_min, Sw_max = 0.57, 0.74

    
    # fig, ax = plt.subplots(2,2,figsize=(10,5))
    # for i in range(len(insts)):
    #     # Clip
    #     Sos_OF_clipped = np.clip(Sos_OF[i], So_min, So_max)
    #     Sgs_OF_clipped = np.clip(Sgs_OF[i], Sg_min, Sg_max)
    #     Sws_OF_clipped = np.clip(Sws_OF[i], Sw_min, Sw_max)

    #     # Rescale to [0,1] for RGB
    #     So_scaled = (Sos_OF_clipped - So_min) / (So_max - So_min)
    #     Sg_scaled = (Sgs_OF_clipped - Sg_min) / (Sg_max - Sg_min)
    #     Sw_scaled = (Sws_OF_clipped - Sw_min) / (Sw_max - Sw_min)

    #     # Build RGB image: So -> red, Sg -> green, Sw -> blue
    #     rgb = np.zeros((ny, nx, 3))
    #     rgb[..., 0] = So_scaled  # red channel
    #     rgb[..., 1] = Sg_scaled  # green channel
    #     rgb[..., 2] = Sw_scaled  # blue channel    

    #     # Plot
    #     ax.flatten()[i].imshow(rgb, origin="lower", aspect="equal")
    #     ax.flatten()[i].set_xlabel("x index")
    #     ax.flatten()[i].set_ylabel("y index")

    # plt.show()

    fig, ax = plt.subplots(2, 2, figsize=(10, 4), constrained_layout=True)  # auto-adjust spacing

    for i in range(len(insts)):
        # # Clip
        # Sos_OF_clipped = np.clip(Sos_OF[i], So_min, So_max)
        # Sgs_OF_clipped = np.clip(Sgs_OF[i], Sg_min, Sg_max)
        # Sws_OF_clipped = np.clip(Sws_OF[i], Sw_min, Sw_max)

        # # Rescale to [0,1] for RGB
        # So_scaled = (Sos_OF_clipped - So_min) / (So_max - So_min)
        # Sg_scaled = (Sgs_OF_clipped - Sg_min) / (Sg_max - Sg_min)
        # Sw_scaled = (Sws_OF_clipped - Sw_min) / (Sw_max - Sw_min)

        # # Build RGB image
        # rgb = np.zeros((ny, nx, 3))
        # rgb[..., 0] = So_scaled  # red
        # rgb[..., 1] = Sg_scaled  # green
        # rgb[..., 2] = Sw_scaled  # blue    

        # # Plot
        # ax_i = ax.flatten()[i]
        # ax_i.imshow(rgb, origin="lower", aspect="equal")
        # ax_i.set_xlabel("x index")
        # ax_i.set_ylabel("y index")

        # # Optional: reduce tick label frequency to save space
        # ax_i.xaxis.set_major_locator(plt.MaxNLocator(4))
        # ax_i.yaxis.set_major_locator(plt.MaxNLocator(4))


        # Suppose you have x and y coordinate arrays
        x = np.linspace(0, 360, 25)  # e.g., x_min=0, x_max=10
        y = np.linspace(0, 100, 5)  # e.g., y_min=0, y_max=5

        # Clip + rescale as before
        Sos_OF_clipped = np.clip(Sos_OF[i], So_min, So_max)
        Sgs_OF_clipped = np.clip(Sgs_OF[i], Sg_min, Sg_max)
        Sws_OF_clipped = np.clip(Sws_OF[i], Sw_min, Sw_max)

        So_scaled = (Sos_OF_clipped - So_min) / (So_max - So_min)
        Sg_scaled = (Sgs_OF_clipped - Sg_min) / (Sg_max - Sg_min)
        Sw_scaled = (Sws_OF_clipped - Sw_min) / (Sw_max - Sw_min)

        rgb = np.zeros((ny, nx, 3))
        rgb[..., 0] = So_scaled
        rgb[..., 1] = Sg_scaled
        rgb[..., 2] = Sw_scaled

        # Plot with extent to map coordinates
        ax_i = ax.flatten()[i]
        ax_i.imshow(rgb, origin="lower", aspect="equal",
                    extent=[x.min(), x.max(), y.min(), y.max()])
        ax_i.set_xlabel("x [m]")
        ax_i.set_ylabel("y [m]")

        # Optional: set ticks on coordinate grid
        ax_i.set_xticks(np.linspace(x.min(), x.max(), 5))  # 5 ticks along x
        ax_i.set_yticks(np.linspace(y.min(), y.max(), 4))  # 4 ticks along y

    plt.savefig(path+'_field.pdf',dpi=300)
    plt.show()












    


