#include "noWell.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    noWell,
    0
);

addToRunTimeSelectionTable
(
    wellModel,
    noWell,
    dictionary
);

noWell::noWell(const dictionary& dict)
:
    wellModel(dict)
{}

void noWell::correct
(
    volScalarField& qa,
        volScalarField& qb,
        const volScalarField& Fa,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod
) const
{

}

} // End namespace Foam
