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
        # return 1.0 - np.exp(-0.5 * sigma**2 - sigma) - DP

    sigma_initial = 0.5
    sigma_solution = fsolve(func, sigma_initial)[0]

    return sigma_solution


def permeability_distribution(DP, k_mean, n_samples=10000):
    """
    Generate permeability distribution from DP coefficient.
    """

    sigma = sigma_from_DP(DP)

    # compute mu from mean permeability
    mu = np.log(k_mean) - 0.5 * sigma**2

    # generate permeability field
    k = np.random.lognormal(mean=mu, sigma=sigma, size=n_samples)

    return k, mu, sigma

# eval k distribution
np.random.seed(1)

DP = 0.8
k_mean = 100.0  # mD
k, mu, sigma = permeability_distribution(DP, k_mean)

plt.hist(k, bins=100)
plt.xlabel("permeability [mD]")
plt.ylabel("frequency")
plt.title(f"Permeability histogram | DP = {DP}")
plt.show()

# generate K

nx = 110
ny = 30

# Elementwise
N = nx*ny
perm = np.random.lognormal(mu, sigma, N)
perm *= 9.869233e-16 # convert mD to m2

# # Layers
# perm_layers = np.random.lognormal(mu, sigma, ny)
# perm_layers *= 9.869233e-16 # convert mD to m2
# perm = np.zeros(nx*ny)
# for i in range(ny):
#     first_in_line = nx*i
#     perm[first_in_line:first_in_line+nx] = perm_layers[i]

plt.hist(perm, bins=100)
plt.xlabel("permeability [m2]")
plt.ylabel("frequency")
plt.title(f"Permeability histogram | DP = {DP}")
plt.show()


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