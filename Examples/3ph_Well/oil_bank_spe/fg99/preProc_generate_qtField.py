L = 2000              # [m]
H = 100               # [m]
W = 750               # [m]
phi = 0.07            # [-] porosity
PV = L * W * H * phi  # [m3] 

Qyear = 0.2 * PV       # [m3/year] 
year = 365 * 24 * 3600  # [s]
Qsec = Qyear / year     # [m3/s]

nx = 110
ny = 30
nz = 1

dx = L / nx        # [m]
dy = H / ny        # [m]
dz = W / nz        # [m]
dv = dx * dy * dz  # [m3]
qsec = Qsec / ny   # [m3/s] per well
qt = qsec / dv     # [m3/s] / [m3] = [1/s]

inj_idx = [i*nx for i in range(ny)]
prod_idx = [i*nx + nx - 1 for i in range(ny)]

n = nx*ny

with open("0/qt", "w") as f:
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
    object      qt;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
""")
    
    # Internal field values
    f.write(f"{n}\n(\n")
    for i in range(n):
        if i in inj_idx:
            f.write(f"    {qt}\n")
        elif i in prod_idx:
            f.write(f"    {-qt}\n")
        else:
            f.write("    0\n")
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
