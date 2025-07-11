import os
import subprocess
from jinja2 import Environment, FileSystemLoader, StrictUndefined
from tqdm import tqdm
import numpy as np

def create_next_state(old_dir, new_dir, ut, deltaT=None):
    """
    Create the next state of a sim by copying the last state 
    from a previous simulation to a new case directory.

    Parameters:
        old_dir (str): Path to the last simulation directory.
        new_dir (str): Path to the new simulation directory.
        ut (tuple): Tuple with the velocity field components (ua, ub, uc).
    """
    
    ua, ub, uc = ut

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

        if uc != None:
            template = env.get_template('Uc')
            with open(f'{new_dir}/0/Uc', 'w') as f:
                f.write(template.render(uc=uc))
        
        if deltaT != None:
            template = env.get_template('controlDict')
            with open(f'{new_dir}/system/controlDict', 'w') as f:
                f.write(template.render(deltaT=deltaT))
    except Exception as e:
        print("Error updating the field:", e)
    # ======================================================================================



if __name__ == '__main__':

    # ======================================================================================
    # Run the first simulation as a starting point, initial conditions should be already set
    first_folder = "fg31.4"
    result = subprocess.run(
        ["bash", "run_solver.sh", first_folder, "> log.txt"],
        check=True, capture_output=True, text=True
    )

    fgs = np.array([
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

    values_ua = np.array([
                3.3180749568552257e-06,
                3.6021914669223394e-06,
                3.8761609587727705e-06,
                4.200865541706615e-06,
                4.505276088207095e-06,
                4.799539616490891e-06,
                5.103950162991371e-06,
                5.428654745925216e-06,
                5.733065292425695e-06,
                6.0882109300095885e-06,
                6.402768494726749e-06,
                6.707179041227228e-06,
                8.462613192713328e-06,
                8.827905848513902e-06,
                9.558491160115054e-06,
                9.609226251198465e-06,
                9.670108360498562e-06,
                9.74113748801534e-06,
                9.791872579098752e-06,
                9.883195743048898e-06,
                9.913636797698944e-06,
                1.0035401016299136e-05,
                1.0065842070949183e-05  
        ])
    
    values_ub = np.array([
                5.81221203451582e-06,
                5.5849188264621285e-06,
                5.3657432329817834e-06,
                5.105979566634708e-06,
                4.862451129434324e-06,
                4.627040306807287e-06,
                4.3835118696069036e-06,
                4.123748203259828e-06,
                3.880219766059444e-06,
                3.5961032559923294e-06,
                3.3444572042186007e-06,
                3.1009287670182177e-06,
                1.6965814458293381e-06,
                1.4043473211888783e-06,
                8.198790719079574e-07,
                7.792909990412283e-07,
                7.305853116011508e-07,
                6.737620095877288e-07,
                6.331739367209985e-07,
                5.601154055608822e-07,
                5.357625618408448e-07,
                4.38351186960691e-07,
                4.139983432406536e-07
        ])
    
    values_uc = np.array([
                1.453053008628955e-06,
                1.3962297066155321e-06,
                1.3414358082454459e-06,
                1.276494891658677e-06,
                1.215612782358581e-06,
                1.1567600767018217e-06,
                1.0958779674017259e-06,
                1.030937050814957e-06,
                9.70054941514861e-07,
                8.990258139980824e-07,
                8.361143010546502e-07,
                7.752321917545544e-07,
                4.2414536145733453e-07,
                3.5108683029721956e-07,
                2.0496976797698934e-07,
                1.9482274976030708e-07,
                1.826463279002877e-07,
                1.684405023969322e-07,
                1.5829348418024963e-07,
                1.4002885139022054e-07,
                1.339406404602112e-07,
                1.0958779674017276e-07,
                1.034995858101634e-07
        ])
    
    deltaTs = 0.1 * np.ones_like(fgs)

    fgs_pairs = np.array([ (fgs[i], fgs[i+1]) for i in range(len(fgs)-1) ])
    print(fgs_pairs)

    for (fg_pair, ua, ub, uc, deltaT) in tqdm(
        zip(fgs_pairs, values_ua[1:], values_ub[1:], values_uc[1:], deltaTs), total=len(fgs_pairs)
    ):  

        new_folder = f"fg{round(float(fg_pair[1]*100),1)}"
        old_folder = f"fg{round(float(fg_pair[0]*100),1)}"

        utt = (ua, ub, uc)

        print(fg_pair, old_folder, new_folder)
        
        create_next_state(old_folder, new_folder, utt, deltaT)

        result = subprocess.run(
            ["bash", "run_solver.sh", new_folder, "> log.txt"],
            check=True, capture_output=True, text=True
        )
    # ======================================================================================
