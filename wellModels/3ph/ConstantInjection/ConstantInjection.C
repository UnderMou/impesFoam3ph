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
    Fa_inj_(readScalar(dict.lookup("Fa_inj"))),
    Fb_inj_(readScalar(dict.lookup("Fb_inj")))
{}

void ConstantInjection::correct
(
    volScalarField& qa,
        volScalarField& qb,
        const volScalarField& Fa,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod
) const
{

    qa = Fa_inj_ * qt_inj - Fa * qt_prod;
    qb = Fb_inj_ * qt_inj - Fb * qt_prod; 

}

} // End namespace Foam
