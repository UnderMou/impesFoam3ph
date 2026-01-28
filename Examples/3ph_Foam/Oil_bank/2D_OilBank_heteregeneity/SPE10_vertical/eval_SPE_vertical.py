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
    # KT.shape[2]//factor, factor
    KT.shape[2]
)
# Take mean over the 2x2 blocks (axes 1 and 3)
# kt = kt.mean(axis=(1, 3, 5))  # Result shape: (110, 30, 85)
kt = kt.mean(axis=(1, 3))  # Result shape: (110, 30, 85)

print(kt.shape)

# vertical layer
ilayer = 15
perm = kt[ilayer, :, ::-1].T.flatten()
# perm = kt[:, ::2, ilayer].flatten()
# perm = kt[:, 1::2, ilayer].flatten()

print(kt[ilayer, :, ::-1].shape)
print(perm.shape)

A = kt[ilayer, :, ::-1].T

plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(figsize=(13, 4))

im = ax.imshow(np.log10(A), aspect='auto', origin='lower', cmap='coolwarm')

# choose how many ticks you want
scale = 10
nx = int(A.shape[1]/(2*scale))
nz = int(A.shape[0]/scale)

# tick positions (indices)
ax.set_xticks(np.linspace(0, A.shape[1]-1, nx))
ax.set_yticks(np.linspace(0, A.shape[0]-1, nz))

# tick labels (physical coordinates)
x = np.linspace(0,300,nx)
z = np.linspace(0,300,nz)
ax.set_xticklabels(np.round(np.linspace(x[0], x[-1], nx), 1))
ax.set_yticklabels(np.round(np.linspace(z[0], z[-1], nz), 1))

ax.set_xlabel(r"$x$ [m]")
ax.set_ylabel(r"$z$ [m]")

plt.colorbar(im, ax=ax, label=r"$\log_{10}(k)$")
plt.tight_layout()
plt.savefig('SPEverticalLayer.pdf', dpi=300)
plt.show()


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
    outlet
    {
        type            zeroGradient;
    }

    inlet
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
}
""")