void correct_Fdry
(
    volScalarField& Fdry,
    const volScalarField Sb, 
    const scalar fmmob,
    const scalar SF,
    const scalar sfbet
)
{    
    const double pi = 3.141592653589793; 
    Fdry = 0.5 + (1.0 / pi) * Foam::atan(sfbet * (Sb - SF));
}

void correct_Nca
(
    volScalarField& Nca,
    const volVectorField U,
    const dimensionedScalar mu_b,
    const scalar sigma_ba
)
{
    Nca = (mu_b * mag(U)) / sigma_ba;
}

void correct_Fshear
(
    volScalarField& Fshear,
    volScalarField& Nca,
    const scalar fmcap,
    const scalar epcap
)
{ 

    forAll(Nca,i){
        if (Nca[i] >= fmcap){
            Fshear[i] = Foam::pow(fmcap/Nca[i], epcap);
        }
    } 

    // Correct at the boundary
    forAll(Nca.boundaryField(), patchi)
    {
        fvPatchScalarField& NcaPatch = Nca.boundaryFieldRef()[patchi];
        fvPatchScalarField& FshearPatch = Fshear.boundaryFieldRef()[patchi];

        forAll(NcaPatch, facei)
        {
            if (NcaPatch[facei] >= fmcap){
                FshearPatch[facei] = Foam::pow(fmcap/NcaPatch[facei], epcap);
            }
        }
    }

}

void correct_MRF
(
    volScalarField& MRF,
    const scalar fmmob,
    const volScalarField Fdry,
    const volScalarField Fshear
)
{
    MRF = 1.0 / (1.0 + fmmob*Fdry*Fshear);
}

void correct_kraMRF
(
    volScalarField& kra,
    const volScalarField MRF
)
{
    kra *= MRF;
}