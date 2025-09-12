scalar Gd = 0.0;
scalar Epsilon = 0.0;
scalar meanGd = 0.0;
scalar meanEpsilon = 0.0;

if (mesh.nInternalFaces())
{
    scalarField Gd_values
        (   
            K*(rho_c.value()-rho_a.value())*mag(g) / (mag(U)*mu_c+VSMALL)
        );

    Gd = gMax(Gd_values);

    meanGd = 0.5*(gSum(Gd_values)/gSum(mesh.V().field()));

    // scalarField Pe_values
    //     (  
    //         mag(U)*mu_c.value()* 
    //         K*(rho_c.value()-rho_a.value())*mag(g) / (mag(U)*mu_c+VSMALL)
    //     );

    // scalarField Epsilon_values
    //     (   
    //         K*(rho_c.value()-rho_a.value())*mag(g) / (mag(U)*mu_c+VSMALL)
    //     );

    // Gd = gMax(Gd_values);

    // meanGd = 0.5*(gSum(Gd_values)/gSum(mesh.V().field()));
}

Info<< "Gd mean: " << meanGd << "," << " max: " << Gd << endl;
// Info<< "Epsilon mean: " << meanEpsilon << "," << " max: " << Epsilon << endl;

// CFLUse = CoNum;
// maxDeltaTFact = maxCo/(CFLUse + SMALL);

// ************************************************************************* //
