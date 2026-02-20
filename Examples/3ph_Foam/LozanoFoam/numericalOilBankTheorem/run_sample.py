import uqtopus as uqt
import numpy as np
import csv
import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from tqdm import tqdm

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
    exp_config = {
        'solver': 'Allrun', 
        'input_path': "templates/S3C2",  
        'output_path': "experiments/init_inj", 
        'parameter_ranges': {},
        'nthreads': 1
    }

    with open("initial_conditions.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        data = [list(map(float, row)) for row in reader]
    X_init = np.array(data)
    X_init = np.delete(X_init, 1, axis=1)
    # print(X_init)

    with open("injection_conditions.csv", "r") as f:
        reader = csv.reader(f)
        next(reader) 
        data = [list(map(float, row)) for row in reader]
    X_inj = np.array(data)
    # X_inj = X_inj[np.argmin(np.linalg.norm(np.array([0.3125, 0, 0.6875]) - X_inj[:,None], axis=2), axis=0)]
    # X_inj = X_inj[np.argmin(np.linalg.norm(np.array([0.05, 0, 0.95]) - X_inj[:,None], axis=2), axis=0)]
    # X_inj = X_inj[np.argmin(np.linalg.norm(np.array([0, 0, 1.0]) - X_inj[:,None], axis=2), axis=0)]
    X_inj = np.delete(X_inj, 1, axis=1)
    # print(X_inj)

    Ni = X_init.shape[0]
    Nj = X_inj.shape[0]
    X = np.hstack([
        np.repeat(X_init, Nj, axis=0),
        np.tile(X_inj, (Ni,1))
    ])
    with open("samples.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sw", "so", "sg"])  # header
        writer.writerows(X)

    # for i in range(len(X)):
    #     exp_config['output_path'] = f"experiments/const_inj/run_{i}"
    #     uqt.run_simulation(params={'0__Sb__Sbinit' : X[i, 0],
    #                             '0__Sa__Sainit' : X[i, 1],
    #                             '0__Sb__Sbinj' : X[i, 2],
    #                             '0__Sa__Sainj' : X[i, 3]},
    #                     exp_config=exp_config, verbose=False)


    n_workers = 5

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(run_case, i, X[i], exp_config)
            for i in range(len(X))
        ]

        for _ in tqdm(as_completed(futures), total=len(futures)):
            pass
        
    