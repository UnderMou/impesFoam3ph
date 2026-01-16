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
    volScalarField& Sc,
    const surfaceScalarField& phia,
    const volScalarField& eps
) const
{

}

} // End namespace Foam
