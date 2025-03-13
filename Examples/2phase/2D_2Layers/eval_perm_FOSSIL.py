import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

# ts = np.linspace(0,150000,601, endpoint=True, dtype=int)
# print(ts.tolist())

nx = 220
ny = 60
perm = np.zeros(nx*ny)
nHalf = int(nx*ny/2)

perm[nHalf:] = 0.1
perm[:nHalf] = 1.0

for i in range(ny):
    first_in_line = nx*i

    perm[first_in_line] = 0.1

with open("perm_2Layers.pbt", "w") as f:
    for value in perm:
        f.write(f"{value}\n")
print("File written!")