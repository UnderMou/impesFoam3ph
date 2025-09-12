scalar CoNum = 0.0;
scalar meanCoNuma = 0.0;
scalar meanCoNumb = 0.0;
scalar meanCoNumc = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhia
        (
            fvc::surfaceSum(mag(phia))().internalField()
        );
    scalarField sumPhib
        (
            fvc::surfaceSum(mag(phib))().internalField()
        );
    scalarField sumPhic
        (
            fvc::surfaceSum(mag(phic))().internalField()
        );

    CoNum = max(0.5*gMax(sumPhia/mesh.V().field())*runTime.deltaTValue(),
                max(
                    0.5*gMax(sumPhib/mesh.V().field())*runTime.deltaTValue(),
                    0.5*gMax(sumPhic/mesh.V().field())*runTime.deltaTValue()
                )
            );

    meanCoNuma = 
        0.5*(gSum(sumPhia)/gSum(mesh.V().field()))*runTime.deltaTValue();
    meanCoNumb = 
        0.5*(gSum(sumPhib)/gSum(mesh.V().field()))*runTime.deltaTValue();
    meanCoNumc = 
        0.5*(gSum(sumPhic)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNuma << "," << meanCoNumb << "," << meanCoNumc << " max: " << CoNum << endl;

// CFLUse = CoNum;
// maxDeltaTFact = maxCo/(CFLUse + SMALL);

// ************************************************************************* //
