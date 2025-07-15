volScalarField correct_pc(const volScalarField Sb){

    // parameters
    scalar Spc_irr = 0;
    scalar Spc_max = 1;
    scalar pc0 = 1e3;
    double alpha = 0.5;

    volScalarField SbSafe = Foam::max(Sb, Spc_irr + SMALL);
    volScalarField Se = (SbSafe-Spc_irr)/(Spc_max-Spc_irr);
    // volScalarField pc = pc0 * Foam::pow(Se, -alpha); 

    volScalarField SeSafe = Foam::max(Se, SMALL);
    volScalarField pc = pc0 * Foam::pow(SeSafe, -alpha); 

    return pc;
}

volScalarField correct_dpcdS(const volScalarField Sb){

    // parameters
    scalar Spc_irr = 0;
    scalar Spc_max = 1;
    scalar pc0 = 1e3;
    double alpha = 0.5;

    volScalarField SbSafe = Foam::max(Sb, Spc_irr + SMALL);
    volScalarField Se = (SbSafe-Spc_irr)/(Spc_max-Spc_irr);
    volScalarField SeSafe = Foam::max(Se, SMALL);


    volScalarField dpcdS = -alpha * pc0 * Foam::pow(SeSafe, -alpha - 1.0) / (Spc_max - Spc_irr);
    return dpcdS;
}