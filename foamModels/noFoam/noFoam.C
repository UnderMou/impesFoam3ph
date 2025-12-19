#include "noFoam.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    noFoam,
    0
);

addToRunTimeSelectionTable
(
    foamModel,
    noFoam,
    dictionary
);

noFoam::noFoam(const dictionary& dict, foamAuxFields* aux)
:
    foamModel(dict, aux)
{}

void noFoam::correct
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

    // if (!aux_ || !aux_->Nca)
    // {
    //     FatalErrorInFunction
    //         << "noFoam requires auxiliary fields but none provided"
    //         << exit(FatalError);
    // }

    // volScalarField& Nca    = *aux_->Nca;
    // volScalarField& Fdry   = *aux_->Fdry;
    // volScalarField& Fshear = *aux_->Fshear;
    // volScalarField& MRF    = *aux_->MRF;

    // correct_Nca(Nca, U); 
    // correct_Fdry(Fdry, Sb);
    // correct_Fshear(Fshear, Nca);
    // correct_MRF(MRF, Fdry, Fshear);
    // kra *= MRF;
}

} // End namespace Foam
