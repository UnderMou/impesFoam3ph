import uqtopus as uqt
import numpy as np
import pandas as pd
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
            'constant__transportProperties__fmmob': X_row.iloc[0],
            'constant__transportProperties__SF': X_row.iloc[1],
            'constant__transportProperties__sfbet': X_row.iloc[2],
            'constant__transportProperties__epcap': X_row.iloc[3],
            'constant__transportProperties__fmoil': X_row.iloc[4],
            'constant__transportProperties__floil': X_row.iloc[5],
            'constant__transportProperties__epoil': X_row.iloc[6],
        },
        exp_config=exp_config,
        verbose=False
    )

if __name__ == "__main__":

    # User defined ===========
    exp_name = "CTscanSA"
    # ========================

    pathCasePrep = f"casesPrep/{exp_name}"
    
    with open("exps_setup.json", "r") as f:
        exp_configs = json.load(f)

    exp_config = exp_configs[exp_name]

    # READING INPUTS
    X = pd.read_csv(f"{pathCasePrep}/input_samples.csv")
    print(X.head())
    print(X.shape)
    # print(X.iloc[0,:])
    # exit()

    # RUN SAMPLES

    n_workers = exp_config['nthreads']

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(run_case, i, X.iloc[i], exp_config)
            for i in range(len(X))
        ]

        for _ in tqdm(as_completed(futures), total=len(futures)):
            pass
        
    
