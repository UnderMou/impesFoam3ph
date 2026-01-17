#include "ConstantInjection.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    ConstantInjection,
    0
);

addToRunTimeSelectionTable
(
    wellModel,
    ConstantInjection,
    dictionary
);

ConstantInjection::ConstantInjection(const dictionary& dict)
:
    wellModel(dict),
    Fb_inj_(readScalar(dict.lookup("Fb_inj")))
{}

void ConstantInjection::correct
(
    volScalarField& qb,
    const volScalarField& Fb, 
    const volScalarField& qt_inj,
    const volScalarField& qt_prod
) const
{

    qb = Fb_inj_ * qt_inj - Fb * qt_prod; 

}

} // End namespace Foam
