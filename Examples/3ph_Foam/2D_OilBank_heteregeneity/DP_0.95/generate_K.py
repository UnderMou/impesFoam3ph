import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def sigma_from_DP(DP):
    """
    Compute lognormal standard deviation (sigma)
    from a given DP coefficient.
    """

    def func(sigma):
        return 1.0 - np.exp(-1.5 * sigma**2 - sigma) - DP

    sigma_initial = 0.5
    sigma_solution = fsolve(func, sigma_initial)[0]

    return sigma_solution


def permeability_from_DP(DP, k_mean, n=10000):
    """
    Generate permeability distribution using
    DP definition from log-normal probability plot.
    """

    # standard deviation in log10 space
    sigma_logk = -np.log10(1.0 - DP)

    # mean in log space
    mu_logk = np.log10(k_mean)

    # sample
    logk = np.random.normal(mu_logk, sigma_logk, n)
    k = 10**logk

    return k, mu_logk, sigma_logk


np.random.seed(1)

DP = 0.95
k_mean = 100.0  # mD
print(f"mean k = {k_mean*9.869233e-16} m2")

k, mu, sigma = permeability_from_DP(DP, k_mean)

print(f"sigma_log10 = {sigma:.3f}")
print(f"mean permeability = {k.mean():.1f} mD")

# Plot (paper-style)
plt.hist(np.log10(k), bins=100)
plt.xlabel("log-permeability [mD]")
plt.ylabel("frequency")
plt.title(f"Permeability histogram | DP = {DP}")
plt.show()

nx = 100
ny = 50

# Elementwise
N = nx*ny
log_perm = np.random.normal(mu, sigma, N)
perm = 10**log_perm
perm *= 9.869233e-16 # convert mD to m2

# # Layers
# log_perm = np.random.normal(mu, sigma, ny)
# perm_layers = 10**log_perm
# perm_layers *= 9.869233e-16 # convert mD to m2
# perm = np.zeros(nx*ny)
# for i in range(ny):
#     first_in_line = nx*i

#     perm[first_in_line:first_in_line+nx] = perm_layers[i]




# WRITE OPENFOAM-LIKE FILE

# Path to write the OpenFOAM-like file
file_path = '0/K'

# Open file and write in the OpenFOAM format
with open(file_path, 'w') as f:
    f.write(
        '''/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  9
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      K;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
'''
    )
    f.write(f"{len(perm)}\n(\n")
    for value in perm:
        f.write(f"{value:.6e}\n")  # Write each perm value in scientific notation
    f.write(")\n")

    f.write(
        ''';
boundaryField
{
    inlet
    {
        type            zeroGradient;
    }
    outlet
    {
        type            zeroGradient;
    }
    walls
    {
        type            zeroGradient;
    }
    emptyWalls
    {
        type            empty;
    }
}'''
    )

print(f"File written to {file_path}")