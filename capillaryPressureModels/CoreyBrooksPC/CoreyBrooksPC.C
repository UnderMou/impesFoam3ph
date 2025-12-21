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
    Spc_irr_(readScalar(dict.lookup("Spc_irr"))),
    Spc_max_(readScalar(dict.lookup("Spc_max"))),
    pc0_(readScalar(dict.lookup("pc0"))),
    alpha_(readScalar(dict.lookup("alpha")))
{

}

void CoreyBrooksPC::correct
(
    volScalarField& pc,
    volScalarField& dpcds,
    const volScalarField& Sb
) const
{   
    // correct dpcds
    volScalarField SbSafe = Foam::max(Sb, Spc_irr_ + SMALL);
    volScalarField Se = (SbSafe-Spc_irr_)/(Spc_max_-Spc_irr_);
    volScalarField SeSafe = Foam::max(Se, SMALL);
    dpcds = -alpha_ * pc0_ * Foam::pow(SeSafe, -alpha_ - 1.0) / (Spc_max_ - Spc_irr_);

    // correct pc
    pc = pc0_ * Foam::pow(SeSafe, -alpha_);
}

} // End namespace Foam
