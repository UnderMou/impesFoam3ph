#include "MehrabiFoam.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    MehrabiFoam,
    0
);

addToRunTimeSelectionTable
(
    foamModel,
    MehrabiFoam,
    dictionary
);

MehrabiFoam::MehrabiFoam(const dictionary& dict, foamAuxFields* aux)
:
    foamModel(dict, aux),
    MRF_(readScalar(dict.lookup("MRF"))),
    Cr_(readScalar(dict.lookup("Cr")))
{}

void MehrabiFoam::correct
(
    volScalarField& kra,
    const volVectorField U,
    const volScalarField& Sa,
    const volScalarField& Sb,
    volScalarField& Sc,
    const surfaceScalarField& phia,
    const volScalarField& eps,
    const volScalarField& K,
    const volVectorField& gradP
) const
{
    // Info << "fmmob = " << fmmob_ << endl;
    // Info << "mu_w = " << mu_b_.value() << endl;

    if (!aux_ || !aux_->MRF)
    {
        FatalErrorInFunction
            << "MehrabiFoam requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    // volScalarField& Nca    = *aux_->Nca;
    // volScalarField& Fdry   = *aux_->Fdry;
    // volScalarField& Fshear = *aux_->Fshear;
    // volScalarField& Fsurf  = *aux_->Fsurf;
    // volScalarField& Foil   = *aux_->Foil;

    volScalarField& MRF    = *aux_->MRF;
    volScalarField& Cs     = *aux_->Cs;

    MRF = 1.0 + (MRF_ - 1.0)*pos(Cs - Cr_); 

    kra /= MRF_;
}

} // End namespace Foam
