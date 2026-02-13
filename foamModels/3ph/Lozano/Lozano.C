#include "Lozano.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    Lozano,
    0
);

addToRunTimeSelectionTable
(
    foamModel,
    Lozano,
    dictionary
);

Lozano::Lozano(const dictionary& dict, foamAuxFields* aux)
:
    foamModel(dict, aux),
    MRF_(readScalar(dict.lookup("MRF")))
{}

void Lozano::correct
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

    // if (!aux_ || !aux_->Fdry)
    // {
    //     FatalErrorInFunction
    //         << "Lozano requires auxiliary fields but none provided"
    //         << exit(FatalError);
    // }

    // volScalarField& Nca    = *aux_->Nca;
    // volScalarField& Fdry   = *aux_->Fdry;
    // volScalarField& Fshear = *aux_->Fshear;
    // volScalarField& Fsurf  = *aux_->Fsurf;
    // volScalarField& Foil   = *aux_->Foil;
    // volScalarField& MRF    = *aux_->MRF;
    // volScalarField& Cs     = *aux_->Cs;

    kra /= MRF_;
}

} // End namespace Foam
