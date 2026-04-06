scalar meanCoNuma = 0.0;
scalar meanCoNumb = 0.0;
scalar maxCoNuma = 0.0;
scalar maxCoNumb = 0.0;

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

    maxCoNuma = 0.5*gMax(sumPhia/mesh.V().field())*runTime.deltaTValue();
    maxCoNumb = 0.5*gMax(sumPhib/mesh.V().field())*runTime.deltaTValue();

    meanCoNuma = 0.5*(gSum(sumPhia)/gSum(mesh.V().field()))*runTime.deltaTValue();
    meanCoNumb = 0.5*(gSum(sumPhib)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Courant Number phase a. mean:" << meanCoNuma << ", max: " << maxCoNuma << 
                      " | phase b. mean:" << meanCoNumb << ", max: " << maxCoNumb << endl;
