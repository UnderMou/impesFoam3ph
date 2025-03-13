import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
plt.rcParams.update({'font.size': 16})

# Openfoam data
data_openfoam = pd.read_csv('2Spots_pressure_openfoam.csv')
of_p = data_openfoam['p'].to_numpy()
of_p /= 1e8
of_x = data_openfoam['arc_length'].to_numpy()
L = of_x[-1]
of_x /= L


# Analytical data
data_analytical = pd.read_csv('2Spots_pressure_carolinaX.csv')
cx_p = data_analytical.iloc[:,1].to_numpy()
cx_x = data_analytical.iloc[:,0].to_numpy()
L = cx_x[-1]
cx_x /= L

plt.figure(figsize=(7,7))
plt.plot(of_x, of_p, label='OpenFoam')
plt.scatter(cx_x, cx_p, label='Analytical')

plt.xlabel(r'$d/D$')
plt.ylabel(r'$p$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("Pressure_comparison.pdf", dpi=300)
plt.show()