import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

# ts = np.linspace(0,150000,601, endpoint=True, dtype=int)
# print(ts.tolist())

nx = 220
ny = 60
perm = np.zeros(nx*ny)
nHalf = int(nx/2)


# counts = int((nx*ny)/(nx/2))


for i in range(ny):
    first_in_line = nx*i

    perm[first_in_line:first_in_line+nHalf] = 0.1
    perm[first_in_line+nHalf:first_in_line+nx] = 1.0

with open("perm_2LayersSeq.pbt", "w") as f:
    for value in perm:
        f.write(f"{value}\n")
print("File written!")