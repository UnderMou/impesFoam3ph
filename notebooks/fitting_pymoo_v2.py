import numpy as np
import os
from pymoo.core.problem import Problem
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.optimize import minimize
from pymoo.termination.default import DefaultSingleObjectiveTermination

from jinja2 import Environment, FileSystemLoader
import subprocess
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd
from aux_funcs import *
from multiprocessing import Pool, Manager

from scipy.linalg import norm

def run_simulation(base_dir, experiment_name, 
                    a_exp = None,
                    b_exp = None,
                    c_exp = None,
                    kra_max = None,
                    krb_max = None, 
                    krc_max = None,
                   verbose=True):
    """
    Runs an OpenFOAM simulation with the given parameters.

    Parameters:
        base_dir (str): Path to the base OpenFOAM case directory.
        experiment_name (str): Name of the experiment to create.
    """

    new_dir = "../simulator/experiments/" + experiment_name

    try:
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        else:
            if verbose:
                print(" -- The directory already exists. Files will be overwritten. --")
            
        result = subprocess.run(
            ["cp", "-a", f'{base_dir}/.', new_dir],
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error copying the files:", e.stderr)

    env = Environment(
        loader=FileSystemLoader(new_dir),
        trim_blocks=True,
        lstrip_blocks=True
    )

    template = env.get_template('constant/transportProperties') # Location of the template file with foam parameters

    output = template.render(a_exp = a_exp,
                             b_exp = b_exp,
                             c_exp = c_exp,
                             kra_max = kra_max,
                             krb_max = krb_max,
                             krc_max = krc_max
                            )

    # Overwrite the template transportProperties file with the new values
    with open(os.path.join(new_dir, 'constant', 'transportProperties'), 'w') as f:
        f.write(output)

    
    result = subprocess.run(
        ["bash", "run_solver.sh", new_dir],  
        check=True,
        capture_output=True,
        text=True
    )
    print(result.stdout)
       

def process_simulation(params, base_case_directory, experiment_name):
    j, (a_exp, b_exp, c_exp, kra_max, krb_max, krc_max) = params
    experiment_name = experiment_name + f"/sample_{j:03d}"
    
    try:
        run_simulation(base_case_directory, experiment_name, 
                        a_exp = a_exp,
                        b_exp = b_exp,
                        c_exp = c_exp,
                        kra_max = kra_max,
                        krb_max = krb_max,
                        krc_max = krc_max,
                        verbose=False)
        # Eval errors
        desired_times = [800]
        ti = 0
        tf = 2000
        write_interval = 200
        t = np.linspace(ti,tf,int(tf/write_interval + 1))
        idt = [np.argmin(np.abs(t-t_find)) for t_find in desired_times]

        # CT Scan data
        tD = [0.36]
        CTscan_results = []
        for i in range(len(tD)):
            Sg_CT_file = 'Tang/Sg_PVI_'+str(tD[i])+'.csv'
            Sw_CT_file = 'Tang/Sw_PVI_'+str(tD[i])+'.csv'
            So_CT_file = 'Tang/So_PVI_'+str(tD[i])+'.csv'
            CTscan_results.append(get_CTscanData(So_CT_file, Sg_CT_file, Sw_CT_file))

        swc = 0.197
        sgr = 0.013
        sor = 0.103

        nGet = -2
        res = []
        for i in range(len(idt)):
            data_dict = parse_openfoam_case('../simulator/experiments/'+experiment_name, variables=['Sa','Sb'])
            Sa = data_dict['Sa'].iloc[idt[i]][:nGet]
            Sb = data_dict['Sb'].iloc[idt[i]][:nGet]
            Sc = 1.0 - Sa - Sb

            # Sa = S_norm(Sa,sgr,swc,sgr,sor)
            # Sb = S_norm(Sb,swc,swc,sgr,sor)
            # Sc = S_norm(Sc,sor,swc,sgr,sor)

            x = np.linspace(0,1,len(Sa))
            L = x[-1]
            x = x/L
            idx = np.argmin(np.abs(x-0.85))

            SoCT = np.interp(x,CTscan_results[i]['x_SoCT'],CTscan_results[i]['So_CT'])
            SgCT = np.interp(x,CTscan_results[i]['x_SgCT'],CTscan_results[i]['Sg_CT'])
            SwCT = np.interp(x,CTscan_results[i]['x_SwCT'],CTscan_results[i]['Sw_CT'])

            error_Sw = np.abs(SwCT - Sb)
            error_Sw = error_Sw[:idx]
            # error_Sg = np.abs(SgCT - Sa)/np.max(SgCT)
            error_So = np.abs(SoCT - Sc)
            error_So = error_So[:idx]

            fig, axs = plt.subplots(1, 2, figsize=(12, 4)) 
            axs[0].plot(x,SoCT,'b')
            axs[0].scatter(CTscan_results[i]['x_SoCT'],CTscan_results[i]['So_CT'],c='k')
            axs[0].plot(x,Sc,'gray')
            axs[0].scatter(x[:idx],error_So,c='r')

            axs[1].plot(x,SwCT,'b')
            axs[1].scatter(CTscan_results[i]['x_SwCT'],CTscan_results[i]['Sw_CT'],c='k')
            axs[1].plot(x,Sb,'gray')
            axs[1].scatter(x[:idx],error_Sw,c='r')
            plt.savefig(f'plots/plot_{j}.png',dpi=200)
            plt.close()

    except subprocess.CalledProcessError as e:
        print("Error running the script:", e.stderr) 

        # Eval errors
        desired_times = [800]
        ti = 0
        tf = 2000
        write_interval = 200
        t = np.linspace(ti,tf,int(tf/write_interval + 1))
        idt = [np.argmin(np.abs(t-t_find)) for t_find in desired_times]

        # CT Scan data
        tD = [0.36]
        CTscan_results = []
        for i in range(len(tD)):
            Sg_CT_file = 'Tang/Sg_PVI_'+str(tD[i])+'.csv'
            Sw_CT_file = 'Tang/Sw_PVI_'+str(tD[i])+'.csv'
            So_CT_file = 'Tang/So_PVI_'+str(tD[i])+'.csv'
            CTscan_results.append(get_CTscanData(So_CT_file, Sg_CT_file, Sw_CT_file))

        swc = 0.197
        sgr = 0.013
        sor = 0.103

        nGet = -2
        res = []
        for i in range(len(idt)):
            data_dict = parse_openfoam_case('../simulator/experiments/'+experiment_name, variables=['Sa','Sb'])
            Sa = data_dict['Sa'].iloc[idt[i]][:nGet]
            Sb = data_dict['Sb'].iloc[idt[i]][:nGet]
            Sc = 1.0 - Sa - Sb

            # Sa = S_norm(Sa,sgr,swc,sgr,sor)
            # Sb = S_norm(Sb,swc,swc,sgr,sor)
            # Sc = S_norm(Sc,sor,swc,sgr,sor)

            x = np.linspace(0,1,len(Sa))
            L = x[-1]
            x = x/L
            idx = np.argmin(np.abs(x-0.85))

            SoCT = np.interp(x,CTscan_results[i]['x_SoCT'],CTscan_results[i]['So_CT'])
            SgCT = np.interp(x,CTscan_results[i]['x_SgCT'],CTscan_results[i]['Sg_CT'])
            SwCT = np.interp(x,CTscan_results[i]['x_SwCT'],CTscan_results[i]['Sw_CT'])

           
            

            error_Sw = 1.0 * np.ones_like(x)
            error_Sw = error_Sw[:idx]
            # error_Sg = np.abs(SgCT - Sa)/np.max(SgCT)
            error_So = 1.0 * np.ones_like(x)
            error_So = error_So[:idx]

            fig, axs = plt.subplots(1, 2, figsize=(12, 4)) 
            axs[0].plot(x,SoCT,'b')
            axs[0].scatter(CTscan_results[i]['x_SoCT'],CTscan_results[i]['So_CT'],c='k')
            axs[0].plot(x,Sc,'gray')
            axs[0].scatter(x[:idx],error_So,c='r')

            axs[1].plot(x,SwCT,'b')
            axs[1].scatter(CTscan_results[i]['x_SwCT'],CTscan_results[i]['Sw_CT'],c='k')
            axs[1].plot(x,Sb,'gray')
            axs[1].scatter(x[:idx],error_Sw,c='r')
            plt.savefig(f'plots/plot_{j}.png',dpi=200)
            plt.close()

    return j, norm(error_So, ord=2) + norm(error_Sw, ord=2)

class CTscan(Problem):

    def __init__(self, **kwargs):

        experiment_name = 'fitting_pymoo'
        BaseCase_dir = '/home/anderson/OpenFOAM/anderson-9/run/uqsa/impesFoam3ph/simulator/base_cases/Lyu_CT_exp_pymoo/'
        new_dir = "../simulator/experiments/" + experiment_name

        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        else:
            print(" -- The directory already exists. Files will be overwritten. --")

        # Parameters ranges
        a_expR      = [1,10] # [1,5] # 1.62 
        b_expR      = [1,10] # [1,5] # 3.86 
        c_expR      = [1,10] # [1,5] # 2.54 
        kra_maxR    = [0.01,1.5] # [0.1,1] # 0.83 
        krb_maxR    = [0.01,1.5] # [0.1,1] # 0.247 
        krc_maxR    = [0.01,1.5] # [0.1,1] # 0.584 
        xl = np.array([a_expR[0] ,b_expR[0] ,c_expR[0] ,kra_maxR[0] ,krb_maxR[0] ,krc_maxR[0]])
        xu = np.array([a_expR[1] ,b_expR[1] ,c_expR[1] ,kra_maxR[1] ,krb_maxR[1] ,krc_maxR[1]])

        self.process_func = partial(
            process_simulation, 
            base_case_directory=BaseCase_dir,
            experiment_name=experiment_name
        )

        super().__init__(n_var=len(xl), n_obj=1, n_ieq_constr=0, xl=xl, xu=xu, **kwargs)

    def _evaluate(self, x, out, *args, **kwargs):
        
        nthreads = 6

        with Manager() as manager:
            output = manager.list([None] * len(x))

            params = list(enumerate(x))

            with Pool(nthreads) as pool:
                for i,result in tqdm(
                    pool.imap_unordered(self.process_func, params),
                    total=len(params), 
                    desc='Running simulations',
                    mininterval=1.0     # Updates at most once per second
                ):
                    output[i] = result

            print(list(output))
            id_min = np.argmin(output)
            print('id_min: ', id_min, ' | f_obj: ', output[id_min])
            print('X : ', x[id_min])
        
            out["F"] = np.array(list(output))

if __name__ == "__main__":
    
    problem = CTscan()

    algorithm = DE(
        pop_size=30,
        sampling= FloatRandomSampling(),
        variant="DE/rand/1/bin",
        CR=0.6,
        dither="vector",
        jitter=False
    )

    termination = DefaultSingleObjectiveTermination(
                        n_max_gen=200,
                        # n_max_evals=400000,
                        period=30)

    ind_run = 10

    for k in range(ind_run):
        seed = k+1

        res = minimize(problem,
                    algorithm,
                    termination,
                    seed=seed,
                    verbose=True)
    
        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))

        file_path = 'pymoo_res/'
        with open(file_path + f'res_seed_{seed}.txt', "w") as f:
            f.write("Best solution found: \n")
            f.write("X = %s\n" % res.X)
            f.write("F = %s\n" % res.F)

