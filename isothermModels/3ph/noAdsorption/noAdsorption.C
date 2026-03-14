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
            << "Isotherm model requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& Cs              = *aux_->Cs;
    volScalarField& dCsEqdCs        = *aux_->dCsEqdCs;
    volScalarField& CsEq            = *aux_->CsEq;

    dCsEqdCs = scalar(0.0) * dCsEqdCs; 
    CsEq = scalar(0.0) * Cs;
}

} // End namespace Foam
