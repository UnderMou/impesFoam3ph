#include "noCapillaryPressure.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    noCapillaryPressure,
    0
);

addToRunTimeSelectionTable
(
    capillaryPressureModel,
    noCapillaryPressure,
    dictionary
);

noCapillaryPressure::noCapillaryPressure(const dictionary& dict)
:
    capillaryPressureModel(dict)
{

}

void noCapillaryPressure::correct
(
    volScalarField& pcwo,
    volScalarField& pcgo,
    volScalarField& dpcwo_dsw,
    volScalarField& dpcwo_dsg,
    volScalarField& dpcgo_dsw,
    volScalarField& dpcgo_dsg,
    const volScalarField& Sa,
    const volScalarField& Sb,
    const volScalarField& zeroField
) const
{   
    // // correct dpcds
    // volScalarField SbSafe = Foam::max(Sb, Spc_irr_ + SMALL);
    // volScalarField Se = (SbSafe-Spc_irr_)/(Spc_max_-Spc_irr_);
    // volScalarField SeSafe = Foam::max(Se, SMALL);
    // dpcds = -alpha_ * pc0_ * Foam::pow(SeSafe, -alpha_ - 1.0) / (Spc_max_ - Spc_irr_);

    // // correct pc
    // pc = pc0_ * Foam::pow(SeSafe, -alpha_);
}

} // End namespace Foam
