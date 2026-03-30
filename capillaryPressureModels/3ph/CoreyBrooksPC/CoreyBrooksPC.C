#include "CoreyBrooksPC.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    CoreyBrooksPC,
    0
);

addToRunTimeSelectionTable
(
    capillaryPressureModel,
    CoreyBrooksPC,
    dictionary
);

CoreyBrooksPC::CoreyBrooksPC(const dictionary& dict)
:
    capillaryPressureModel(dict),
    epsilon_(readScalar(dict.lookup("epsilon")))
{

}

void CoreyBrooksPC::correct
(
    volScalarField& pcow,
    volScalarField& pcgo,
    volScalarField& dpcow_dsw,
    volScalarField& dpcow_dsg,
    volScalarField& dpcgo_dsw,
    volScalarField& dpcgo_dsg,
    const volScalarField& Sa,
    const volScalarField& Sb,
    const volScalarField& zeroField
) const
{   
    // correct pcwo
    pcow = -5.0 * epsilon_ * (scalar(2.0) - Sb) * (scalar(1.0) - Sb);
    dpcow_dsw = -5.0 * epsilon_ * (scalar(2.0)*Sb - scalar(3.0));
    dpcow_dsg = zeroField;

    // correct pcgo
    pcgo = 5.0 * epsilon_ * (scalar(2.0) - Sa) * (scalar(1.0) - Sa);
    dpcgo_dsg = 5.0 * epsilon_ * (scalar(2.0)*Sa - scalar(3.0));
    dpcgo_dsw = zeroField;
}

} // End namespace Foam
