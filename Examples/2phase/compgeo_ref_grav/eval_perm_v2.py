import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ts = np.linspace(0,150000,601, endpoint=True, dtype=int)
# print(ts.tolist())
nx = 220
ny = 60
 
perm = np.loadtxt('/home/lamap/OpenFOAM/lamap-9/run/impesFoam3ph-multisolver/Examples/2phase/spe10_layer_36.pbt')
d2m2 = 9.869233e-13
# perm = np.log10(d2m2 * perm)
perm = d2m2 * perm
perm = perm.reshape((60,220))

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
    f.write(f"{nx*ny}\n(\n")
    for row in perm:
        for value in row:
            f.write(f"{value}\n")  # Write each perm value in scientific notation
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
