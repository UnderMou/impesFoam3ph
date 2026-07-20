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
    plt.rcParams.update({'font.size': 22})

    ti = 0
    tf = 4e8
    write_interval = 5e6
    t = np.linspace(ti, tf, int(tf/write_interval)+1)

    nx, ny = 110, 30
    lx, ly = 360, 100

    dx = lx/nx
    dy = ly/ny
    dz = 10.0
    dv = dx*dy*dz
    qt = 4.883307965499746e-08
    Qt = ny * qt * dv

    phi = 0.07
    PV = lx * ly * 10.0 * phi

    PVI = Qt/PV * t

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)

    dirs = ['fg30/', 'fg70/', 'fg99/']

    fgs = [0.3, 0.7, 1.0]

    qo = np.zeros((len(fgs),len(t)))
    qg = np.zeros_like(qo)

    for i, d in enumerate(dirs):

        data_dict = parse_openfoam_case(d, variables=['qa','qb','qt'])
        data_dict.index = (data_dict.index.astype(float) / write_interval).astype(int)

        qas_OF = data_dict['qa']
        qbs_OF = data_dict['qb']
        qts_OF = data_dict['qt']

        insts_id = np.array([np.argmin(np.abs(t-inst)) for inst in t], dtype=int)
        insts_id = insts_id[1:]

        qa_OF = [np.reshape(qas_OF[j], (ny, nx)) for j in insts_id]
        qb_OF = [np.reshape(qbs_OF[j], (ny, nx)) for j in insts_id]
        qt_OF = [np.reshape(qts_OF[j], (ny, nx)) for j in insts_id]
        qo_OF = [q - qa - qb for qa, qb, q in zip(qa_OF, qb_OF, qt_OF)]
        
        for j in insts_id[:-1]:
            qo[i,j] = np.sum(qo_OF[j][:,-1]*dv) # integrate over cell volume and sum ny cells
            qg[i,j] = np.sum(qa_OF[j][:,-1]*dv)
        # axes[0].plot(t[:-1],-qo[i,:-1],linewidth=2,label=r"$f_g^{\textrm{inj.}} = $"+f" {str(fgs[i]*100)}\%")
        axes[0].plot(PVI[:-1],-qg[i,:-1]/-qo[i,:-1],linewidth=2,label=r"$f_g^{\textrm{inj.}} = $"+f" {str(fgs[i])}")
    
    # axes[0].set_ylabel(f"Oil-cut [PV/sec.]")
    axes[0].set_ylabel(f"GOR [$q_g/q_o$]")
    # axes[0].set_xlabel(f"time [sec.]")
    axes[0].set_xlabel(f"injected pore volume [PVI]")
    axes[0].grid()
    axes[0].legend()


    dt = t[1:] - t[:-1]
    Qo = np.zeros_like(qo)
    OOIP = 0.46
    for i in range(len(fgs)):
        for j in insts_id[:-1]:
            Qo[i,j] = np.sum(-qo[i,:j]*dt[:j])
        
        # axes[1].plot(t[:-1],((Qo[i,:-1]/PV)/OOIP)*100,linewidth=2,label=r"$f_g^{\textrm{inj.}} = $"+f" {str(fgs[i]*100)}\%")
        axes[1].plot(PVI[:-1],((Qo[i,:-1]/PV)/OOIP)*100,linewidth=2,label=r"$f_g^{\textrm{inj.}} = $"+f" {str(fgs[i]*100)}\%")
    
    axes[1].set_ylabel(f"ORF [\%]")
    # axes[1].set_xlabel(f"time [sec.]")
    axes[1].set_xlabel(f"injected pore volume [PVI]")
    axes[1].grid()



    plt.savefig('fgs_oilCutORF.pdf', dpi=300, bbox_inches='tight')
    plt.show()
        
        