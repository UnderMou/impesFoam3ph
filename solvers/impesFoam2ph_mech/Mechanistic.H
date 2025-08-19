volScalarField correct_nD_LE(const volScalarField Sb,
                     const scalar Sw_ast,
                     const scalar A
                    )
{   
    
    volScalarField nD_LE = Sb;    // just for initialization
    forAll(nD_LE,i){
        if (Sb[i] <= Sw_ast){
            nD_LE[i] = 0.0;
        }
        else{
            nD_LE[i] = Foam::tanh(A * (Sb[i] - Sw_ast));
        }
    }

    return nD_LE;
}

volScalarField correct_rgrc(const volScalarField Sb,
                            const volScalarField nD,
                            const scalar Kc,
                            const scalar nMAX,
                            const scalar Sw_ast,
                            const scalar A
                           )
{
    volScalarField nD_LE = correct_nD_LE(Sb,Sw_ast,A);
    volScalarField rgrc = Kc * nMAX * (nD_LE - nD);

    return rgrc;
}

volScalarField correct_MRF(const volScalarField nD,
                           const scalar Cmrf
                          )
{
    volScalarField MRF = Cmrf * nD + scalar(1.0);
    return MRF;
}