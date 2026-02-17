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

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams['font.size'] = 16

    # # Reading Eduardo Abreu data
    # EduAbreu_Sg_Grav = np.array([[0, 0.0128, 0.0304, 0.0488, 0.06720000000000001, 0.09040000000000001, 0.1184, 0.15760000000000002, 0.19440000000000002, 0.2544, 0.2888, 0.3168, 0.3432, 0.35760000000000003, 0.3672, 0.37520000000000003, 0.3808, 0.38480000000000003, 0.388, 0.39440000000000003, 0.40240000000000004, 0.4112, 0.4192, 0.4248, 0.42800000000000005, 0.42960000000000004, 0.4304, 0.432, 0.43360000000000004, 0.43760000000000004, 0.44480000000000003, 0.4632, 0.48000000000000004, 0.4968, 0.5008, 0.5032, 0.504, 0.5064000000000001, 0.5088, 0.5112, 0.516, 0.536, 0.5992000000000001, 0.664, 0.7608, 0.8816, 0.9984000000000001],
    #                              [0.2729032258064516, 0.2696774193548387, 0.26129032258064516, 0.2535483870967742, 0.24838709677419354, 0.23935483870967741, 0.23096774193548386, 0.22064516129032258, 0.21096774193548387, 0.1967741935483871, 0.18774193548387097, 0.18064516129032257, 0.17225806451612904, 0.16709677419354838, 0.16193548387096773, 0.15548387096774194, 0.14709677419354839, 0.13806451612903226, 0.13161290322580646, 0.12580645161290321, 0.12580645161290321, 0.12645161290322582, 0.12838709677419355, 0.13225806451612904, 0.1464516129032258, 0.16709677419354838, 0.1858064516129032, 0.20129032258064516, 0.2103225806451613, 0.2135483870967742, 0.21483870967741936, 0.21483870967741936, 0.21548387096774194, 0.21548387096774194, 0.2103225806451613, 0.2, 0.18903225806451612, 0.17419354838709677, 0.16322580645161291, 0.15548387096774194, 0.15096774193548387, 0.1503225806451613, 0.1503225806451613, 0.1503225806451613, 0.1503225806451613, 0.1503225806451613, 0.1503225806451613]])
    # # EduAbreu_Sw_Grav = pd.read_csv('EduAbreu_Sw_noGrav.csv')
    # EduAbreu_So_Grav = np.array([[-0.002352941176470588, 0.01764705882352941, 0.02941176470588235, 0.05058823529411764, 0.08352941176470588, 0.1176470588235294, 0.15294117647058822, 0.18941176470588234, 0.2364705882352941, 0.28823529411764703, 0.31999999999999995, 0.3423529411764706, 0.36, 0.3717647058823529, 0.3776470588235294, 0.38352941176470584, 0.3905882352941176, 0.40352941176470586, 0.40352941176470586, 0.4152941176470588, 0.42117647058823526, 0.42470588235294116, 0.4270588235294117, 0.4294117647058823, 0.4305882352941176, 0.4317647058823529, 0.4317647058823529, 0.4329411764705882, 0.4341176470588235, 0.4341176470588235, 0.43882352941176467, 0.4529411764705882, 0.47764705882352937, 0.4952941176470588, 0.5023529411764706, 0.5047058823529411, 0.5070588235294117, 0.508235294117647, 0.5129411764705882, 0.5247058823529411, 0.5729411764705882, 0.6388235294117647, 0.7070588235294117, 0.7847058823529411, 0.8729411764705881, 0.9941176470588234],
    #                              [0.016209773539928488, 0.040047675804529205, 0.0629320619785459, 0.09439809296781883, 0.13539928486293207, 0.17163289630512515, 0.20405244338498213, 0.23456495828367105, 0.27079856972586414, 0.3079856972586412, 0.33087008343265795, 0.3489868891537545, 0.36519666269368295, 0.37854588796185934, 0.39380214541120384, 0.41096543504171634, 0.4224076281287247, 0.42526817640047676, 0.42526817640047676, 0.42431466030989273, 0.4300357568533969, 0.4443384982121573, 0.4901072705601907, 0.5549463647199047, 0.5988081048867699, 0.6388557806912992, 0.6912991656734208, 0.7094159713945173, 0.7218116805721096, 0.7284862932061978, 0.734207389749702, 0.734207389749702, 0.735160905840286, 0.735160905840286, 0.7418355184743742, 0.7589988081048867, 0.7790226460071514, 0.7904648390941598, 0.798092967818832, 0.8, 0.8, 0.8, 0.8009535160905841, 0.8, 0.8, 0.8009535160905841]])
 
    # Define t vector
    ti = 0
    tf = 1
    write_interval = 0.001
    t = np.linspace(ti,tf,int(tf/write_interval)+1)
    print(t)
    t[-1] = int(1)
    # time_dirs = [0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06,0.0625,0.065,0.0675,0.07,0.0725,0.075,0.0775,0.08,0.0825,0.085,0.0875,0.09,0.0925,0.095,0.0975,0.1,0.1025,0.105,0.1075,0.11,0.1125,0.115,0.1175,0.12,0.1225,0.125]
    time_dirs = [f'{x:.3f}'.rstrip('0').rstrip('.') if x % 1 else f'{int(x)}' for x in t]
    print(time_dirs)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    

    # Reading OpenFOAM files
    path = ['CaseB', 'CaseB', 'CaseB_noG', 'CaseB_noG']
    colors = ['b','b','k', 'k']
    styles = ['-','-','--','--']
    t_find = [0.275, 0.495, 0.275, 0.495]

    for i in range(len(path)):
        data_dict = parse_openfoam_case(path[i], variables=['Sb','Sa'], time_dirs=time_dirs)
        if i < 2: EduAbreu_Sg_Grav = pd.read_csv('Sg_t'+str(i+1)+'.csv')
        if i < 2: EduAbreu_So_Grav = pd.read_csv('So_t'+str(i+1)+'.csv')
        
        id_t = np.argmin(np.abs(t-t_find[i]))
        print(t[id_t])
        Sw = data_dict['Sb'].iloc[id_t]
        Sg = data_dict['Sa'].iloc[id_t]
        So = 1.0 - Sw - Sg

        # x domain 
        x = np.linspace(0,1,len(Sw))
        
        # Plot profiles
        axes[0].plot(x, So, c=colors[i], linewidth=1, linestyle=styles[i], label=r'impesFoam')
        if i < 2: axes[0].plot(EduAbreu_So_Grav.iloc[:,0], EduAbreu_So_Grav.iloc[:,1], linestyle='-.', c='r', label=r'de Abreu, E. (2007)')
        # axes[0].scatter(EduAbreu_So_Grav[0,:], EduAbreu_So_Grav[1,:], s=2, c=colors[2], label=r'de Abreu, E. (2007)')
        axes[0].grid()
        axes[0].set_xlabel(r"$x$")
        axes[0].set_ylabel(r"$S_o$")
        # axes[0].set_ylim([0,0.85])
        # axes[0].set_xlim([0,0.55])

        axes[1].plot(x, Sg, c=colors[i], linewidth=1, linestyle=styles[i], label=r'impesFoam')
        if i < 2: axes[1].plot(EduAbreu_Sg_Grav.iloc[:,0], EduAbreu_Sg_Grav.iloc[:,1], linestyle='-.', c='r', label=r'de Abreu, E. (2014)')
        # axes[1].scatter(EduAbreu_Sg_Grav[0], EduAbreu_Sg_Grav[1], s=2, c=colors[2], label=r'de Abreu, E. (2007)')
        axes[1].grid()
        axes[1].set_xlabel(r"$x$")
        axes[1].set_ylabel(r"$S_g$")
        # axes[1].set_ylim([0.,0.4])
    
    # axes[1].legend(loc='upper right')

    fluid_legend = [
        Line2D([0], [0], color='b', linestyle='-', lw=1, label='ImpesFoam - $G_d=0.0475$'),
        Line2D([0], [0], color='k', linestyle='--', lw=2, label='ImpesFoam - $G_d=0$'),
        Line2D([0], [0], color='r', linestyle='-.', lw=2, label='de Abreu, E. (2014)')
    ]
    legend1 = axes[1].legend(handles=fluid_legend, loc="upper right")
    axes[1].add_artist(legend1)

    plt.tight_layout()
    # fig.suptitle("gravity effects")
    plt.savefig('EduardoAbreuComparison_Gd_0.085.pdf', dpi=300)
    plt.show()