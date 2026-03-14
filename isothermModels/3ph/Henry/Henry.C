#include "Henry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    Henry,
    0
);

addToRunTimeSelectionTable
(
    isothermModel,
    Henry,
    dictionary
);

Henry::Henry(const dictionary& dict, foamAuxFields* aux)
:
    isothermModel(dict, aux),
    Kd_(readScalar(dict.lookup("Kd")))
{}

void Henry::correct
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
    const volScalarField& oneField  = *aux_->oneField;

    dCsEqdCs = Kd_ * oneField;
    CsEq = Kd_ * Cs;

}

} // End namespace Foam
