import numpy as np
import pandas as pd
import os
import re
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

if __name__ == '__main__':

    k = 1.9541081340000003e-12
    ut = 1.058334e-05
    L = 0.4 

    plt.style.use('science')
    plt.rcParams.update({'font.size': 18})

    # muApp_alg = np.array([0.59969355, 0.62642292, 0.64667623, 0.67002855, 0.69647233,
    #    0.71830841, 0.74470437, 0.76795762, 0.7897291 , 0.8164127 ,
    #    0.84063516, 0.86535766, 0.97458892, 0.98712777, 0.76443993, 0.74391864,
    #    0.75324925, 0.713967  , 0.68651833, 0.62874674, 0.5616681, 0.47864273, 0.44107242])
    # fg_alg = np.array([0.32882345, 0.35525579, 0.37584728, 0.40017983, 0.42848233,
    #    0.45243118, 0.4820718 , 0.50879242, 0.53431572, 0.5662656 ,
    #    0.5959106 , 0.62681496, 0.77726838, 0.80075283, 0.90481342, 0.90808511,
    #    0.90660429, 0.91277518, 0.91700633, 0.92577732, 0.93585736, 0.94838858, 0.95414779])
    # sorted_muApp = [mu for fg, mu in sorted(zip(fg_alg, muApp_alg), key=lambda pair: pair[0])]
    # sorted_fg = [mu for fg, mu in sorted(zip(fg_alg, muApp_alg), key=lambda pair: pair[0])]
    # muApp_alg = sorted_muApp
    # fg_alg = sorted_fg
    
    fg_exp = np.array([
            0.31351869606903165,
            0.34036433365292423,
            0.36625119846596355,
            0.3969319271332694,
            0.42569511025886864,
            0.45349952061361454,
            0.48226270373921376,
            0.5129434324065196,
            0.5417066155321189,
            0.5752636625119847,
            0.6049856184084371,
            0.6337488015340363,
            0.7996164908916586,
            0.8341323106423777,
            0.9031639501438159,
            0.9079578139980824,
            0.9137104506232022,
            0.9204218600191754,
            0.9252157238734419,
            0.9338446788111218,
            0.9367209971236816,
            0.9482262703739214,
            0.9511025886864811
        ])
    
    muApp_exp = [
            0.6133682830930537,
            0.6408912188728701,
            0.6566186107470511,
            0.673656618610747,
            0.6998689384010484,
            0.7169069462647444,
            0.7444298820445608,
            0.7614678899082568,
            0.7771952817824377,
            0.7994757536041939,
            0.8217562254259501,
            0.8492791612057666,
            0.8768020969855831,
            0.8217562254259501,
            0.783748361730013,
            0.7391874180865006,
            0.6959370904325033,
            0.6513761467889907,
            0.618610747051114,
            0.5792922673656618,
            0.5360419397116645,
            0.49803407601572736,
            0.45347313237221487
        ]
    PD_exp = (np.array(muApp_exp) * ut * L / k)/1e6

    folders = [f for f in os.listdir('.') if os.path.isdir(f) and re.match(r'fg\d+', f)]
    folders_sorted = sorted(folders, key=lambda x: float(re.search(r'\d+(\.\d+)?', x).group()))

    print(folders_sorted)

    

    Swc = 0.135
    Sgr = 0.2
    Sor = 0.1
    phi = 0.22

    PD = []
    muApp = []
    muApp_IMPES_steady = []


    plt.figure(figsize=(15,4))

    total_t = [0]
    tD = [0]
    for i in range(len(folders_sorted)):
    
        ti = 0
        tf = 40000
        write_interval = 500
        t = np.linspace(ti,tf,int(tf/write_interval + 1))
        nAvoid = -1
        time_steps = len(t[1:])

        path = folders_sorted[i]
        data_dict = parse_openfoam_case(path, variables=['p'])

        last_t = total_t[-1]
        for j in range(1, time_steps):

            p_out = data_dict['p'].iloc[j][:nAvoid][-1]
            p_in = data_dict['p'].iloc[j][:nAvoid][0]

            PD.append(p_in - p_out)
            muApp.append(k*PD[-1]/(ut*L))
            total_t.append(float(t[j] + last_t))
            tD.append(float(t[j] + last_t)*ut/((1.0-Swc-Sgr-Sor)*L*phi))
        print(total_t[-1], muApp_exp[i])

        muApp_IMPES_steady.append(muApp[-1])

        tf_adm = tf * ut / ((1.0-Swc-Sgr-Sor)*L*phi)

        # plt.axvline(x=total_t[-1], color='gray', alpha=0.7, linestyle='--')
        plt.axvline(x=tD[-1], color='gray', alpha=0.7, linestyle='--')

        # plt.text(total_t[-1]-tf*0.9, max(muApp_exp)*1.23, folders_sorted[i][2:] + r'\%', rotation=0, color='k', fontsize=14)
        plt.text(tD[-1]-tf_adm*0.9, max(PD_exp)*1.23, folders_sorted[i][2:] + r'\%', rotation=0, color='k', fontsize=14)

        # plt.text(tD[-1]-tf*0.9, max(muApp_exp)*1.23, folders_sorted[i][2:] + r'\%', rotation=0, color='k', fontsize=14)

        # plt.scatter([total_t[-1]], [muApp_exp[i]], c='k')
        # plt.scatter([total_t[-1]], [muApp_alg[i]], c='r')
    
    # plt.scatter([total_t[-1]], [muApp_exp[-1]], c='k', label = 'Steady-state experimental data')
    # plt.scatter([total_t[-1]], [muApp_alg[-1]], c='r', label = 'Steady-state Algebraic model')
    
    print(total_t)

    
    tD = (ut*np.array(total_t)) / ((1.0-Swc-Sgr-Sor)*L*phi)

    # plt.plot(total_t[:-1], muApp, label = 'ImpesFOAM')
    # plt.plot(total_t[:-1], np.array(PD)/1e6, label = 'ImpesFOAM')
    plt.plot(tD[:-1], np.array(PD)/1e6, label = 'ImpesFOAM')
    plt.grid(axis='y')
    # plt.ylabel(r'$\mu_{app}$ [Pa.s]')
    plt.ylabel(r'$\Delta p$ [MPa]')
    plt.xlabel(r'Dimensionless time [PVI]')
    plt.tight_layout()
    # plt.ylim([-0.05,1.3])
    # plt.legend(loc='best')
    plt.savefig('DP_openFOAM.pdf', dpi=500)
    plt.show()
    plt.close()

    model_alg = pd.read_csv('model_alg.csv')
    df_sorted = model_alg.sort_values(by='fg_alg', ascending=True)
    df_sorted = df_sorted.reset_index(drop=True)

    
    fg_alg = df_sorted['fg_alg']
    id = np.argmin(np.abs(0.1-fg_alg))
    muApp_alg = df_sorted['muApp_alg']

    plt.figure(figsize=(6,6))
    plt.plot(fg_alg[id:], muApp_alg[id:], c='r', label = 'Steady-state algebraic model')
    plt.scatter(fg_exp, muApp_exp, c='k', label = 'Steady-state experimental data')
    plt.scatter(fg_exp, muApp_IMPES_steady, c='b', label = 'Steady-state ImpesFOAM')
    plt.grid(True)
    plt.grid(which='minor', linestyle=':', linewidth=0.5, color='gray')
    plt.ylabel(r'$\mu_{app}$ [Pa.s]')
    plt.xlabel(r'$f_g$ [-]')
    plt.legend(loc='center left', bbox_to_anchor=(0.01, 1.2))
    plt.tight_layout()
    plt.savefig('steadyState_muApp_fg_openFOAM.pdf', dpi=300)
    plt.show()