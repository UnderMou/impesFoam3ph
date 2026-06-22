import uqtopus as uqt
import numpy as np
import csv
import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from tqdm import tqdm
import json

def run_case(i, X_row, exp_config_base):

    exp_config = copy.deepcopy(exp_config_base)
    exp_config['output_path'] += f"/run_{i}"

    uqt.run_simulation(
        params={
            '0__Sb__Sbinit': X_row[0],
            '0__Sa__Sainit': X_row[1],
            '0__Sb__Sbinj':  X_row[2],
            '0__Sa__Sainj':  X_row[3]
        },
        exp_config=exp_config,
        verbose=False
    )

if __name__ == "__main__":

    # User defined ===========
    exp_name = "Lyu_SFL"
    # ========================

    pathCasePrep = f"casesPrep/{exp_name}"
    
    with open("exps_setup.json", "r") as f:
        exp_configs = json.load(f)

    exp_config = exp_configs[exp_name]

    # READING INPUTS
    with open(f"{pathCasePrep}/initial_conditions_{exp_name}.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        data = [list(map(float, row)) for row in reader]
    X_init = np.array(data)
    X_init = np.delete(X_init, 1, axis=1)   # delete So column
    print(X_init[:5])

    with open(f"{pathCasePrep}/injection_conditions_{exp_name}.csv", "r") as f:
        reader = csv.reader(f)
        next(reader) 
        data = [list(map(float, row)) for row in reader]
    X_inj = np.array(data)
    X_inj = np.delete(X_inj, 1, axis=1) # delete So column
    print(X_inj[:5])

    # DEFINING SAMPLES SET FOR SIMULATIONS
    Ni = X_init.shape[0]
    Nj = X_inj.shape[0]
    X = np.hstack([
        np.repeat(X_init, Nj, axis=0),
        np.tile(X_inj, (Ni,1))
    ])
    print(X[:5])

    with open(exp_config['output_path']+"/samples.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sw", "sg", "sw_inj", "sg_inj"])  # header
        writer.writerows(X)
    print(X.shape)
    print(X[:5])

    # RUN SAMPLES

    n_workers = exp_config['nthreads']

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(run_case, i, X[i], exp_config)
            for i in range(len(X))
        ]

        for _ in tqdm(as_completed(futures), total=len(futures)):
            pass
        
    
