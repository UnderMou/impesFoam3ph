// volScalarField correct_pc(const volScalarField Sb){

//     volScalarField pc = Foam::min(Sb, SMALL);

//     return pc;
// }

// volScalarField correct_dpcdS(const volScalarField Sb){

//     volScalarField dpcdS = Foam::min(Sb, SMALL);
//     return dpcdS;
// }

volScalarField correct_pc(const volScalarField Sb){

    // parameters
    scalar Swc = 0.2;
    scalar Sgr = 0.18;
    scalar pc0 = 330.0;
    double alpha = 0.01;

    volScalarField pc = pc0 * Foam::pow(scalar(1.0) - Sb - Sgr, alpha) / (Sb - Swc + VSMALL); 

    return pc;
}

volScalarField correct_dpcdS(const volScalarField Sb){

    // parameters
    scalar Swc = 0.2;
    scalar Sgr = 0.18;
    scalar pc0 = 330.0;
    double alpha = 0.01;

    volScalarField dpcdS = ((pc0 * alpha * (-scalar(1.0) + Swc + Sgr)) / Foam::pow(Sb - Swc + VSMALL, 2))*
                            Foam::pow((scalar(1.0) - Sb - Sgr)/(Sb - Swc),alpha-1); 

    return dpcdS;
}