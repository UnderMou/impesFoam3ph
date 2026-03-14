scalar Gd = 0.0;
scalar meanGd = 0.0;
// scalar Epsilon = 0.0;
// scalar meanEpsilon = 0.0;

if (mesh.nInternalFaces())
{
    // scalarField Gd_values
    //     (   
    //         K*(rho_b.value()-rho_a.value())*mag(g) / (mag(U)*mu_b+VSMALL)
    //     );

    const scalarField Umag = mag(U);
    const scalar Umin = 1e-12;

    scalarField Gd_values
    (
        K*(rho_b.value()-rho_a.value())*mag(g)
        / ( max(mag(U), scalar(Umin)) * mu_b )
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
