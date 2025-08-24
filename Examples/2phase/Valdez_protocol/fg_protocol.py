import os
import subprocess
from jinja2 import Environment, FileSystemLoader, StrictUndefined
from tqdm import tqdm
import numpy as np

def create_next_state(old_dir, new_dir, ut, deltaT=None, endTime=None):
    """
    Create the next state of a sim by copying the last state 
    from a previous simulation to a new case directory.

    Parameters:
        old_dir (str): Path to the last simulation directory.
        new_dir (str): Path to the new simulation directory.
        ut (tuple): Tuple with the velocity field components (ua, ub, uc).
    """
    
    ua, ub = ut

    # Get the last time directory
    last_time = max([int(d) for d in os.listdir(old_dir) if d[0].isdigit()])

    # ======================================================================================
    # Copy the files from the last state to the new directory
    try:
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        else:
            # print(f" -- The directory '{new_dir}' already exists. Files will be overwritten. --")
            subprocess.run(
                ["rm", "-rf", new_dir],
                check=True, capture_output=True, text=True
            )
            os.makedirs(new_dir)

        os.makedirs(f'{new_dir}/constant')
        os.makedirs(f'{new_dir}/system')
        os.makedirs(f'{new_dir}/0')

        subprocess.run(
            ["cp", "-a", f'{old_dir}/constant/.', f'{new_dir}/constant'],
            check=True, capture_output=True, text=True
        )
        subprocess.run(
            ["cp", "-a", f'{old_dir}/system/.', f'{new_dir}/system'],
            check=True, capture_output=True, text=True
        )
        subprocess.run(
            ["cp", "-a", f'{old_dir}/{last_time}/.', f'{new_dir}/0'],
            check=True, capture_output=True, text=True
        )

    except subprocess.CalledProcessError as e:
        print("Error copying the files:", e.stderr)
    # ======================================================================================
    

    # ======================================================================================
    # Update the velocity field based on the contents in 'templates' directory
    env = Environment(loader=FileSystemLoader('templates'), undefined=StrictUndefined)
    
    try:
        if ua != None:
            template = env.get_template('Ua')
            with open(f'{new_dir}/0/Ua', 'w') as f:
                f.write(template.render(ua=ua))
        
        if ub != None:
            template = env.get_template('Ub')
            with open(f'{new_dir}/0/Ub', 'w') as f:
                f.write(template.render(ub=ub))
        
        if deltaT != None:
            template = env.get_template('controlDict')
            with open(f'{new_dir}/system/controlDict', 'w') as f:
                f.write(template.render(deltaT=deltaT, endTime=endTime))
        
    except Exception as e:
        print("Error updating the field:", e)
    # ======================================================================================



if __name__ == '__main__':

    # ======================================================================================
    # Run the first simulation as a starting point, initial conditions should be already set
    first_folder = "fg50.0"
    result = subprocess.run(
        ["bash", "run_solver.sh", first_folder, "> log.txt"],
        check=True, capture_output=True, text=True
    )

    fgs = np.array([
            0.5,
            0.6,
            0.7,
            0.75,
            0.8,
            0.9
        ])
    
    endTimes = np.array([
        10521,
        7238,
        8019,
        10512,
        5095
    ])

    ut = 1.45e-5    # [m/s]

    values_ua = ut * fgs
    
    values_ub = ut - values_ua
    
    deltaTs = 0.1 * np.ones_like(fgs)
    print(deltaTs)

    fgs_pairs = np.array([ (fgs[i], fgs[i+1]) for i in range(len(fgs)-1) ])
    print(fgs_pairs)

    for (fg_pair, ua, ub, deltaT, endTime) in tqdm(
        zip(fgs_pairs, values_ua[1:], values_ub[1:], deltaTs, endTimes), total=len(fgs_pairs)
    ):  

        new_folder = f"fg{round(float(fg_pair[1]*100),1)}"
        old_folder = f"fg{round(float(fg_pair[0]*100),1)}"

        utt = (ua, ub)

        print(fg_pair, old_folder, new_folder)
        
        create_next_state(old_folder, new_folder, utt, deltaT, endTime)

        result = subprocess.run(
            ["bash", "run_solver.sh", new_folder, "> log.txt"],
            check=True, capture_output=True, text=True
        )
    # ======================================================================================
