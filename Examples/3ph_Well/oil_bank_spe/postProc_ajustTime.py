import os
import numpy as np

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# t_find = 8250

# Time controls:
ti = 0
tf = 3e8 
write_interval = 1e6
t = np.linspace(ti,tf,int(tf/write_interval + 1))
# idt = np.argmin(np.abs(t-t_find))

chang_dirs = ['fg30/', 'fg70/', 'fg99/']
base_dir = os.getcwd()

for i in range(len(chang_dirs)):

    direc = os.path.join(base_dir,chang_dirs[i])

    items = os.listdir(direc)
    numerical_folders = [
        item for item in items
        if os.path.isdir(os.path.join(direc, item))
        and is_number(item)
    ]

    numerical_folders = sorted(
        [item for item in items
        if os.path.isdir(os.path.join(direc, item)) and is_number(item)],
        key=lambda x: float(x)
    )

    int_numerical_folders = [int(float(s)) for s in numerical_folders]
    adjustedTime = int_numerical_folders

    print(numerical_folders)
    print(adjustedTime)


    # print(adjustedTime)
    # exit()
    for i in range(len(numerical_folders)):
            # original_num = int(folder)
            # new_num = original_num * scaling_factor + offset
            new_folder_name = str(adjustedTime[i])
            
            original_path = os.path.join(direc, str(numerical_folders[i]))
            new_path = os.path.join(direc, new_folder_name)
            
            # Rename the folder
            os.rename(original_path, new_path)
            print(f'Renamed "{numerical_folders[i]}" to "{new_folder_name}"')

    os.chdir(base_dir)