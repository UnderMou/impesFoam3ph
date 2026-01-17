import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

filename = "spe_perm.dat"

perm_data = np.loadtxt(filename).flatten()
perm_data = perm_data * 9.869233e-13   # D to m2

nx, ny, nz = 60, 220, 85
total_cells = nx * ny * nz

kx = perm_data[:total_cells]
ky = perm_data[total_cells:2*total_cells]
kz = perm_data[2*total_cells:3*total_cells]

# Fortran-like ordering
kx = kx.reshape((nx, ny, nz), order='F')
ky = ky.reshape((nx, ny, nz), order='F')
kz = kz.reshape((nx, ny, nz), order='F')

KT = kx**2 # + ky**2 # + kz*2
KT = np.sqrt(KT)

# DO NOTHING
# kt = KT

# Reshape into blocks of size 2x2 in axes 0 and 1
factor = 2
kt = KT.reshape(
    KT.shape[0]//factor, factor,   # (110, 2)
    KT.shape[1]//factor, factor,   # (30, 2)
    KT.shape[2]          # (85)
)
# Take mean over the 2x2 blocks (axes 1 and 3)
kt = kt.mean(axis=(1, 3))  # Result shape: (110, 30, 85)

print(kt.shape)

ilayer = 35
Delta_layer = 10
# kt_aux = kt[:, ::-1, ilayer-Delta_layer:ilayer+Delta_layer]
# kt_aux = kt[:, ::-1, ilayer:ilayer+Delta_layer]
kt_aux = kt[:, :, ilayer:ilayer+Delta_layer]

for i in range(kt_aux.shape[2]):
    if i == 0: perm = kt_aux[:, ::-1, i].flatten()
    else:
        perm = np.concatenate((perm, kt_aux[:, ::-1, i].flatten()))

print(perm.shape)

perm_df = pd.DataFrame(perm)
perm_df.to_csv(f'spe10_layer{ilayer+1:02}.csv', index=False, header=False)

n = perm.shape[0]

with open("0/K", "w") as f:
    # Header
    f.write("""/*--------------------------------*- C++ -*----------------------------------*\\
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
""")
    
    # Internal field values
    f.write(f"{n}\n(\n")
    for i in range(n):
        f.write(f"    {perm[i]}\n")
    f.write(");\n\n")

    # Boundary field
    f.write("""boundaryField
{
    walls
    {
        type            zeroGradient;
    }

    emptyWalls
    {
        type            empty;
    }
}
""")