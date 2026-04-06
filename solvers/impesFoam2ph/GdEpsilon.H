scalar Gd = 0.0;
scalar meanGd = 0.0;
// scalar Epsilon = 0.0;
// scalar meanEpsilon = 0.0;

if (mesh.nInternalFaces())
{
    scalarField Gd_values
        (   
            K*(rho_b.value()-rho_a.value())*mag(g) / (mag(U)*mu_b+VSMALL)
        );

    Gd = gMax(Gd_values);

    meanGd = average(Gd_values);

    // if capPressModel->type() != "noCapillaryPressure"
    // {
    //     scalarField Pe_values
    //     (  
    //         mag(U)*mu_b.value()* / (pc*K+VSMALL)
    //     );
    // }
    
}

Info<< "Gd mean: " << meanGd << "," << " max: " << Gd << endl;

// ************************************************************************* //
