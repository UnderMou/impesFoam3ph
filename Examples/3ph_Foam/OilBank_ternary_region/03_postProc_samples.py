import uqtopus as uqt
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import ternary
import matplotlib
import pandas as pd
import json

def oilBank_criteria(
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
        return False

    # Locate extrema
    i_min = np.argmin(g)
    i_max = np.argmax(g)

    # Get region between extrema
    i1, i2 = sorted([i_min, i_max])
    Sc_segment = Sc[i1:i2+1]

    # Condition 2: saturation in between must exceed Soi
    if not np.all(Sc_segment > Soi):
        return False

    # Condition 3: saturation in between must exceed Soi, in at least one position
    # if np.max(Sc_segment) < Soi:
       #return False
    
    # Condition 4: saturations at oil bank segment should be a plateau
    # Sc_segment = np.asarray(Sc_segment)
    # tol = 1e-6 * np.max(Sc_segment)
    # if not np.allclose(Sc_segment, Sc_segment[0], atol=tol):
    #     return False

    return True 

def checkCase_OBcondition(Sc, Soi, i, exp_name):
    """
        This function marks the simulation samples that has lead to oscilations 
        in the numerical results. So we ignore them by applying the exclusions
        shown bellow. These exclusion criteria was defined looking at the outputs
        and manually pinned for each case. 
    """

    if exp_name == "CTscan":
        return oilBank_criteria(Sc, Soi) and Soi < 0.6 # not all Soi > 0.6 oscillated, just to make it easier
    
    elif exp_name == "Lozano":
        return oilBank_criteria(Sc, Soi)
    
    elif exp_name == "Lyu":
        return oilBank_criteria(Sc, Soi) and i not in [94, 95, 119, 127, 1071] 
    elif exp_name == "Lyu_fmmob_x2":
        return oilBank_criteria(Sc, Soi) and i not in [95,103,111,119,127,596,598,599,918,919,1071]
    elif exp_name == "Lyu_fmmob_xhalf":
        return oilBank_criteria(Sc, Soi) and i not in [95]
    elif exp_name == "Lyu_fmoil_High":
        return oilBank_criteria(Sc, Soi) and X[i, 0]>0.1
    elif exp_name == "Lyu_fmoil_xhalf":
        return oilBank_criteria(Sc, Soi) and X[i, 0]>0.1 and Soi<=0.2 and i not in [221, 223, 885, 887, 1027, 1291, 1515]
    
    else:
        raise Exception("'checkOil_bank' condition not defined for this 'exp_name'!")


if __name__ == "__main__":

    # User defined ===========
    exp_name = "CTscan_moreInj"
    # ========================

    pathCasePrep = f"casesPrep/{exp_name}"
    
    with open("exps_setup.json", "r") as f:
        exp_configs = json.load(f)

    exp_config = exp_configs[exp_name]

    # Reading the simulation results    =====================================

    n_folders = sum(
        os.path.isdir(os.path.join(exp_config['output_path'], item))
        for item in os.listdir(exp_config['output_path'])
    )
    print(n_folders)

    with open(exp_config['output_path']+"/samples.csv", "r") as f:
        reader = csv.reader(f)
        next(reader) 
        data = [list(map(float, row)) for row in reader]
    X = np.array(data)
    print(X.shape)

    ppt = exp_config['postProcTime']
    start, end, step = ppt['start'], ppt['end'], ppt['step']

    # Check for oil bank formation for all samples  =========================
    fg_inj = []
    checkOil_bank = []
    for i in range(n_folders):
        exp_dir = exp_config['output_path'] + "/run_" + str(i)
        data_single = uqt.parse_openfoam_case(
            case_dir=exp_dir,
            variables=["Sa", "Sb", "Sc", "Fa"],
            time_dirs=[f"{i:g}" for i in np.arange(start, end, step)]   
        )
        Sc = data_single.Sc[-1,:]

        Soi = 1 - X[i, 0] - X[i, 1]
        Sgi = X[i, 1]

        Fa = np.asarray(data_single.Fa[-1,:])
        fg_inj.append(Fa[0])

        # Check for oil bank formation along all simulation time
        temp_check = []
        for j in range(len(data_single.time)):
            Sc = data_single.Sc[j,:]

            if checkCase_OBcondition(Sc, Soi, i, exp_name):

                temp_check.append(True)
                
                print(X[i,0],Soi,Sgi,X[i,2],X[i,3])
                plt.plot(Sc, c='b',alpha=0.5)
                break
            else:
                temp_check.append(False)

        checkOil_bank.append(np.any(temp_check))
        print(i, checkOil_bank[-1])


    # SAVING THE RESULTS

    # 1. Save Oil saturation profiles samples
    plt.savefig(exp_config['output_path']+"/samples_profiles.png", dpi=300)
    plt.show()
    plt.close()

    # 2. Save DataFrame of oil bank formation
    df = pd.DataFrame(X, columns=['sw','sg','sw_inj','sg_inj'])
    df["oil_bank"] = checkOil_bank 
    df["fg_inj"] = fg_inj
    df.to_csv(exp_config['output_path']+"/X.csv", index=False)

    # 3. Save ternary diagram of oil bank formation cases
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
