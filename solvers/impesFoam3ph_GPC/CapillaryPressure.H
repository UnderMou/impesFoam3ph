volScalarField correct_pcwo(const volScalarField Sb, const scalar epsilon){

    volScalarField pcwo = 5.0 * epsilon * (scalar(2.0) - Sb) * (scalar(1.0) - Sb); 

    return pcwo;
}

volScalarField correct_pcgo(const volScalarField Sa, const scalar epsilon){

    volScalarField pcgo = 5.0 * epsilon * (scalar(2.0) - Sa) * (scalar(1.0) - Sa);

    return pcgo;
}

volScalarField correct_Dpcwo_DSw(const volScalarField Sb, const scalar epsilon){

    volScalarField Dpcwo_DSw = 5.0 * epsilon * (scalar(2.0)*Sb - scalar(3.0));
    
    return Dpcwo_DSw;
}

volScalarField correct_Dpcgo_DSg(const volScalarField Sa, const scalar epsilon){

    volScalarField Dpcgo_DSg = 5.0 * epsilon * (scalar(2.0)*Sa - scalar(3.0));
    
    return Dpcgo_DSg;
}