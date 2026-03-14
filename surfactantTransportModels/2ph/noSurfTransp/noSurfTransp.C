#include "noSurfTransp.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    noSurfTransp,
    0
);

addToRunTimeSelectionTable
(
    surfactantTransportModel,
    noSurfTransp,
    dictionary
);

noSurfTransp::noSurfTransp(const dictionary& dict, foamAuxFields* aux)
:
    surfactantTransportModel(dict, aux)
{

}

void noSurfTransp::correct
(
    const volScalarField& Sb,
    const surfaceScalarField& phib,
    const volScalarField& eps
) const
{
    
}

} // End namespace Foam
