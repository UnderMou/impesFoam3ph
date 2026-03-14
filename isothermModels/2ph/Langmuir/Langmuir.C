#include "Langmuir.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    Langmuir,
    0
);

addToRunTimeSelectionTable
(
    isothermModel,
    Langmuir,
    dictionary
);

Langmuir::Langmuir(const dictionary& dict, foamAuxFields* aux)
:
    isothermModel(dict, aux),
    qo_(readScalar(dict.lookup("qo"))),
    KL_(readScalar(dict.lookup("KL")))
{}

void Langmuir::correct
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

    dCsEqdCs = qo_ * KL_ / (Foam::pow(scalar(1.0) + KL_ * Cs,2));
    CsEq = qo_ * KL_ * Cs / (scalar(1.0) + KL_ * Cs);

}

} // End namespace Foam
