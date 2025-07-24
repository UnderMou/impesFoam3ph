
// LINEAR ADSORPTION 

volScalarField correct_CsEq(const volScalarField Cs, const scalar Kd){

    volScalarField CsEq = Kd*Cs;
    
    return CsEq;
}

scalar correct_dCsEqdCs(scalar Kd){

    return Kd;
}
