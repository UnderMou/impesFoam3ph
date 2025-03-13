import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

# DATA SPE10/FOSSIL
spe10 = pd.read_csv('spe10_layer_36.pbt', header=None)
spe10 = spe10.iloc[:,0].to_numpy()
perm_base = 1e-12

# DATA OPENFOAM
openf = pd.read_csv('out.csv', header=None)
openf_x = openf.iloc[:,1].to_numpy()
openf_y = openf.iloc[:,2].to_numpy()
K = perm_base * spe10

nx = 220
ny = 60
for i in range(ny):
    K[i*nx] = 1e-13

# PLOT
# plt.scatter(spe10_x, spe10_y, c=spe10_perm, marker='o', label="Original Data", cmap='viridis', edgecolor='black')
plt.scatter(openf_x, openf_y, c=K, marker='x', label="Interpolated Data", cmap='viridis')
plt.colorbar(label='Permeability')
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Permeability Interpolation')
plt.show()

# WRITE OPENFOAM-LIKE FILE

# Path to write the OpenFOAM-like file
file_path = 'K'

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
    f.write(f"{len(K)}\n(\n")
    for value in K:
        f.write(f"{value:.6e}\n")  # Write each perm value in scientific notation
    f.write(")\n")

    f.write(
        ''';
boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 1e-13;
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