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

STARS::STARS(const dictionary& dict, foamAuxFields* aux)
:
    foamModel(dict, aux),
    fmmob_(readScalar(dict.lookup("fmmob"))),
    SF_(readScalar(dict.lookup("SF"))),
    sfbet_(readScalar(dict.lookup("sfbet"))),
    fmcap_(readScalar(dict.lookup("fmcap"))),
    epcap_(readScalar(dict.lookup("epcap"))),
    sigma_ba_(readScalar(dict.lookup("sigma_ba"))),
    fmsurf_(readScalar(dict.lookup("fmsurf"))),
    epsurf_(readScalar(dict.lookup("epsurf"))),
    FshearActive_(dict.parent().lookupOrDefault<Switch>("FshearActive",false)),
    FsurfActive_(dict.parent().lookupOrDefault<Switch>("FsurfActive",false))
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
    const volVectorField U
) const
{
    Nca = (mu_b_.value() * mag(U)) / sigma_ba_;
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
        } else {
            Fshear[i] = 1.0;
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
            } else {
                FshearPatch[facei] = 1.0;
            }
        }
    }
}

void STARS::correct_Fsurf
(
    volScalarField& Fsurf,
    volScalarField& Cs
) const
{
    // Correct at cell center
    forAll(Cs,i){
        if (Cs[i] < fmsurf_){
            Fsurf[i] = Foam::pow(Cs[i]/fmsurf_, epsurf_);
        } else {
            Fsurf[i] = 1.0;
        }
    } 

    // Correct at the boundaries
    forAll(Cs.boundaryField(), patchi)
    {
        fvPatchScalarField& CsPatch = Cs.boundaryFieldRef()[patchi];
        fvPatchScalarField& FsurfPatch = Fsurf.boundaryFieldRef()[patchi];

        forAll(CsPatch, facei)
        {
            if (CsPatch[facei] < fmsurf_){
                FsurfPatch[facei] = Foam::pow(CsPatch[facei]/fmsurf_, epsurf_);
            } else {
                FsurfPatch[facei] = 1.0;
            }
        }
    }
}

void STARS::correct_MRF
(
    volScalarField& MRF,
    const volScalarField Fdry,
    const volScalarField Fshear,
    const volScalarField Fsurf
) const
{
    MRF = 1.0 / (1.0 + fmmob_*Fdry*Fshear*Fsurf);
}

void STARS::correct
(
    volScalarField& kra,
    const volVectorField U,
    const volScalarField& Sa,
    const volScalarField& Sb,
    const surfaceScalarField& phia,
    const volScalarField& eps
) const
{
    // Info << "fmmob = " << fmmob_ << endl;
    // Info << "mu_w = " << mu_b_.value() << endl;

    if (!aux_ || !aux_->Fdry)
    {
        FatalErrorInFunction
            << "STARS requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& Nca    = *aux_->Nca;
    volScalarField& Fdry   = *aux_->Fdry;
    volScalarField& Fshear = *aux_->Fshear;
    volScalarField& Fsurf  = *aux_->Fsurf;
    volScalarField& MRF    = *aux_->MRF;
    volScalarField& Cs     = *aux_->Cs;

    correct_Fdry(Fdry, Sb);
    if (FshearActive_)
    {
        correct_Nca(Nca, U);
        correct_Fshear(Fshear, Nca);
    }
    if (FsurfActive_)
    {
        correct_Fsurf(Fsurf, Cs);
    }
    correct_MRF(MRF, Fdry, Fshear, Fsurf);
    kra *= MRF;
}

} // End namespace Foam
