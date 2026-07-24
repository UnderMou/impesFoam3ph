nx = 1000

qt =  1e-2  # [1 / s]

idx = 0
first_value = qt
last_value = -qt

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
    f.write(f"{nx}\n(\n")
    for i in range(nx):
        if i == idx:
            f.write(f"    {first_value}\n")
        elif i == nx - 1 - idx:
            f.write(f"    {last_value}\n")
        else:
            f.write("    0\n")
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

    emptyWalls
    {
        type            empty;
    }
}
""")
