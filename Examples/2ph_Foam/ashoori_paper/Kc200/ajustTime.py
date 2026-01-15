import os
import numpy as np

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

t_find = 8250

# Time controls:
ti = 0
tf = 2e6 
write_interval = 1000
t = np.linspace(ti,tf,int(tf/write_interval + 1))
idt = np.argmin(np.abs(t-t_find))


base_dir = '.'
items = os.listdir(base_dir)
numerical_folders = [
    item for item in items
    if os.path.isdir(os.path.join(base_dir, item))
    and is_number(item)
]

numerical_folders = sorted(
    [item for item in items
     if os.path.isdir(os.path.join(base_dir, item)) and is_number(item)],
    key=lambda x: float(x)
)

int_numerical_folders = [int(float(s)) for s in numerical_folders]
adjustedTime = int_numerical_folders

print(numerical_folders)
print(adjustedTime)


for i in range(len(numerical_folders)):
        # original_num = int(folder)
        # new_num = original_num * scaling_factor + offset
        new_folder_name = str(adjustedTime[i])
        
        original_path = os.path.join(base_dir, str(numerical_folders[i]))
        new_path = os.path.join(base_dir, new_folder_name)
        
        # Rename the folder
        os.rename(original_path, new_path)
        print(f'Renamed "{numerical_folders[i]}" to "{new_folder_name}"')


# for i in range(len(numerical_folders)):
#         if str(adjustedTime[i])[-1] == '9':
        
#             new_folder_name = str(adjustedTime[i]+1)
        
#             original_path = os.path.join(base_dir, str(numerical_folders[i]))
#             new_path = os.path.join(base_dir, new_folder_name)
        
#             # Rename the folder
#             os.rename(original_path, new_path)
#             print(f'Renamed "{numerical_folders[i]}" to "{new_folder_name}"')
#         else:
#             print(f'Not renamed "{numerical_folders[i]}"')