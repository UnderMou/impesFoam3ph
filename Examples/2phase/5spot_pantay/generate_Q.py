nx = 51
n = nx**2
qt = 0.0031868

first_value = qt
last_value = -qt

with open("Q", "w") as f:
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
    object      Q;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
""")
    
    # Internal field values
    f.write(f"{n}\n(\n")
    for i in range(n):
        if i == 0:
            f.write(f"    {first_value}\n")
        elif i == n - 1:
            f.write(f"    {last_value}\n")
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
