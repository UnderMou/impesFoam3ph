import uqtopus as uqt
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import ternary
import matplotlib

def check_oilBank(
    Sc,
    Soi,
    grad_min_threshold = -0.01,
    grad_max_threshold =  0.01
):
    Sc = np.asarray(Sc)

    # Compute gradient
    g = np.gradient(Sc)

    g_min = np.min(g)
    g_max = np.max(g)

    # Condition 1: gradient must indicate a bank structure
    if not (g_min < grad_min_threshold and g_max > grad_max_threshold):
        # print("Condition 1 failed")
        return False

    # Locate extrema
    i_min = np.argmin(g)
    i_max = np.argmax(g)

    # Get region between extrema
    i1, i2 = sorted([i_min, i_max])
    Sc_segment = Sc[i1:i2+1]

    # Condition 2: saturation in between must exceed Soi
    if not np.all(Sc_segment > Soi):
        # print("Condition 2 failed")
        return False
    
    # # Condition 3: saturations at oil bank segment should be a plateau
    # Sc_segment = np.asarray(Sc_segment)
    # tol = 1e-6 * np.max(Sc_segment)
    # if not np.allclose(Sc_segment, Sc_segment[0], atol=tol):
    #     print("Condition 3 failed")
    #     return False

    # print("Oil bank detected")
    return True 

exp_config = {
        'solver': 'Allrun', 
        'input_path': "templates/S3C2",  
        'output_path': "experiments/init_inj_stars", 
        'parameter_ranges': {},
        'nthreads': 1
    }

n_folders = sum(
    os.path.isdir(os.path.join(exp_config['output_path'], item))
    for item in os.listdir(exp_config['output_path'])
)
print(n_folders)

with open("samples.csv", "r") as f:
    reader = csv.reader(f)
    next(reader) 
    data = [list(map(float, row)) for row in reader]
X = np.array(data)
print(X)

checkOil_bank = []
for i in range(n_folders):
    exp_dir = exp_config['output_path'] + "/run_" + str(i)
    data_single = uqt.parse_openfoam_case(
        case_dir=exp_dir,
        variables=["Sa", "Sb", "Sc"],
        time_dirs=[f"{i:g}" for i in np.arange(1000, 20000, 1000)]
    )
    Sc = data_single.Sc[-1,:]
    # print(Sc)
    Soi = 1 - X[i, 0] - X[i, 1]
    # print(Soi)
    Sgi = X[i, 1]

    
    # checkOil_bank.append(np.any(Sc > 1.01*Soi))
    # checkOil_bank.append(check_oilBank(Sc, Soi))
    # checkOil_bank.append(np.any((Sc > 1.01*Soi) & (Sgi > 0.5) & (Soi < 0.2)))

    temp_check = []
    for j in range(len(data_single.time)):
        Sc = data_single.Sc[j,:]
        if check_oilBank(Sc, Soi):
            temp_check.append(True)

            # print(j,"|",len(data_single.time))
            plt.plot(Sc, c='b',alpha=0.5)
            break
        else:
            temp_check.append(False)
    checkOil_bank.append(np.any(temp_check))
    print(i, checkOil_bank[-1])

#     if checkOil_bank[-1]:
#         plt.plot(Sc, c='b',alpha=0.5)
plt.show()
plt.close()

# Graph
scale = 1
matplotlib.rcParams['figure.figsize'] = (8, 7)
figure, tax = ternary.figure(scale=scale)
tax.gridlines(multiple=0.1, linewidth=0.5, color="blue")
tax.boundary(linewidth=2.0)
fontsize = 14
tax.right_corner_label("W", fontsize=fontsize, color="blue", fontweight="bold")
tax.top_corner_label("O", fontsize=fontsize, color="black", fontweight="bold")
tax.left_corner_label("G", fontsize=fontsize, color="green", fontweight="bold")

# ALL POINTS
# Initial conditions
Sw, Sg = X[:,0], X[:,1]
So = 1 - Sw - Sg
points = [(sw, so, sg) for sw, so, sg in zip(Sw, So, Sg)]
tax.scatter(points, marker='s', color='dimgray', label="Samples")
# Injection conditions
Sw, Sg = X[:,2], X[:,3]
So = 1 - Sw - Sg
points = [(sw, so, sg) for sw, so, sg in zip(Sw, So, Sg)]
tax.scatter(points, marker='s', color='dimgray')

# OIL BANK
# Initial conditions
Sw, Sg = X[checkOil_bank,0], X[checkOil_bank,1]
So = 1 - Sw - Sg
points = [(sw, so, sg) for sw, so, sg in zip(Sw, So, Sg)]
tax.scatter(points, marker='s', color='red', label="Oil bank")
# Injection conditions
Sw, Sg = X[checkOil_bank,2], X[checkOil_bank,3]
So = 1 - Sw - Sg
points = [(sw, so, sg) for sw, so, sg in zip(Sw, So, Sg)]
tax.scatter(points, marker='s', color='darkred', label="Injection for oil bank")



tax.ticks(axis='lbr', linewidth=1, offset=0.025, multiple=0.1, tick_formats="%.1f")
tax.get_axes().axis('off')
tax.clear_matplotlib_ticks()

plt.legend()

plt.savefig(exp_config['output_path']+"/samples_ternary.png", dpi=300)
# tax.show()