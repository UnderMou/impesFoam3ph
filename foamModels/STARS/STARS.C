#include "STARS.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    STARS,
    0
);

addToRunTimeSelectionTable
(
    foamModel,
    STARS,
    dictionary
);

STARS::STARS(const dictionary& dict)
:
    foamModel(dict),
    fmmob_(readScalar(dict.lookup("fmmob"))),
    SF_(readScalar(dict.lookup("SF"))),
    sfbet_(readScalar(dict.lookup("sfbet"))),
    fmcap_(readScalar(dict.lookup("fmcap"))),
    epcap_(readScalar(dict.lookup("epcap"))),
    sigma_ba_(readScalar(dict.lookup("sigma_ba")))
{}

void STARS::correct_Fdry
(
    volScalarField& Fdry,
    const volScalarField Sb
) const
{
    const double pi = 3.141592653589793; 
    Fdry = 0.5 + (1.0 / pi) * Foam::atan(sfbet_ * (Sb - SF_));
}

void STARS::correct_Nca
(
    volScalarField& Nca,
    const volVectorField U,
    const dimensionedScalar mu_b
) const
{
    Nca = (mu_b.value() * mag(U)) / sigma_ba_;
}

void STARS::correct_Fshear
(
    volScalarField& Fshear,
    volScalarField& Nca
) const
{
    // Correct at cell center
    forAll(Nca,i){
        if (Nca[i] >= fmcap_){
            Fshear[i] = Foam::pow(fmcap_/Nca[i], epcap_);
        }
    } 

    // Correct at the boundaries
    forAll(Nca.boundaryField(), patchi)
    {
        fvPatchScalarField& NcaPatch = Nca.boundaryFieldRef()[patchi];
        fvPatchScalarField& FshearPatch = Fshear.boundaryFieldRef()[patchi];

        forAll(NcaPatch, facei)
        {
            if (NcaPatch[facei] >= fmcap_){
                FshearPatch[facei] = Foam::pow(fmcap_/NcaPatch[facei], epcap_);
            }
        }
    }
}

void STARS::correct_MRF
(
    volScalarField& MRF,
    const volScalarField Fdry,
    const volScalarField Fshear
) const
{
    MRF = 1.0 / (1.0 + fmmob_*Fdry*Fshear);
}

void STARS::correct
(
    volScalarField& kra,
    volScalarField& MRF,
    volScalarField& Fdry,
    volScalarField& Fshear,
    volScalarField& Nca,
    const volVectorField U,
    const volScalarField& Sb
) const
{
    // Info << "fmmob = " << fmmob_ << endl;
    // Info << "mu_w = " << mu_b_.value() << endl;

    correct_Nca(Nca, U, mu_b_); 
    correct_Fdry(Fdry, Sb);
    correct_Fshear(Fshear, Nca);
    correct_MRF(MRF, Fdry, Fshear);
    kra *= MRF;
}

} // End namespace Foam
