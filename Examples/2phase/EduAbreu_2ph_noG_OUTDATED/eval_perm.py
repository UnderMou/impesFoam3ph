import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

# DATA SPE10/FOSSIL
spe10 = pd.read_csv('perm.csv')
spe10_perm = spe10['Permeability'].to_numpy()
spe10_perm = np.power(10, spe10_perm)
# spe10_perm = np.exp(spe10_perm)
spe10_x = spe10['Points:0'].to_numpy()
spe10_y = spe10['Points:1'].to_numpy()

# DATA OPENFOAM
openf = pd.read_csv('out.csv', header=None)
openf_x = openf.iloc[:,2].to_numpy()
openf_y = openf.iloc[:,3].to_numpy()
openf_perm = openf.iloc[:,1].to_numpy()

# # INTERPOLATE
# points = np.column_stack((spe10_x, spe10_y))
# points_new = np.column_stack((openf_x, openf_y))
# perm_new = interp.griddata(points, spe10_perm, points_new, method='linear')
# # print(len(perm_new))

# 2 LAYERS
perm_new = np.zeros_like(openf_perm)
mask_UL = openf_x > 4.5
mask_LL = openf_x <= 4.5
perm_new[mask_UL] = 1.0e-12
perm_new[mask_LL] = 1.0e-13

# # PLOT
# plt.scatter(spe10_x, spe10_y, c=spe10_perm, marker='o', label="Original Data", cmap='viridis', edgecolor='black')
# plt.scatter(openf_x, openf_y, c=perm_new, marker='x', label="Interpolated Data", cmap='viridis')
# plt.colorbar(label='Permeability')
# plt.legend()
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('2D Permeability Interpolation')
# plt.show()


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
    f.write(f"{len(perm_new)}\n(\n")
    for value in perm_new:
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