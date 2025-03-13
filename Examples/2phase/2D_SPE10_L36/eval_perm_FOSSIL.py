import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

perm = pd.read_csv('spe10_layer_36.csv', header=None)
# print(perm.head())
perm = perm.iloc[:,0].to_numpy()
# print(perm)
# print(np.max(perm),np.min(perm))
nx = 220
ny = 60

for i in range(ny):
    perm[i*nx] = 0.1


with open("perm.pbt", "w") as f:
    for value in perm:
        f.write(f"{value}\n")
print("File written!")