#include "noAdsorption.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    noAdsorption,
    0
);

addToRunTimeSelectionTable
(
    isothermModel,
    noAdsorption,
    dictionary
);

noAdsorption::noAdsorption(const dictionary& dict, foamAuxFields* aux)
:
    isothermModel(dict, aux)
{}

void noAdsorption::correct
(
    
) const
{

    if (!aux_ || !aux_->dCsEqdCs)
    {
        FatalErrorInFunction
            << "Isother model requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& Cs              = *aux_->Cs;
    volScalarField& dCsEqdCs        = *aux_->dCsEqdCs;
    volScalarField& CsEq            = *aux_->CsEq;

    dCsEqdCs = 0.0 * dCsEqdCs; // because dCsEqdCs is initially
                               // filled of 1.0
    CsEq = 0.0 * Cs;

}

} // End namespace Foam
