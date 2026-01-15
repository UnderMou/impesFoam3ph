#include "ashooriPC.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    ashooriPC,
    0
);

addToRunTimeSelectionTable
(
    capillaryPressureModel,
    ashooriPC,
    dictionary
);

ashooriPC::ashooriPC(const dictionary& dict)
:
    capillaryPressureModel(dict),
    Swc_(readScalar(dict.lookup("Swc"))),
    Sgr_(readScalar(dict.lookup("Sgr"))),
    pc0_(readScalar(dict.lookup("pc0"))),
    c_(readScalar(dict.lookup("c")))
{

}

void ashooriPC::correct
(
    volScalarField& pc,
    volScalarField& dpcds,
    const volScalarField& Sb
) const
{   
    // correct dpcds
    // volScalarField SbSafe = Foam::max(Sb, SMALL);
    // volScalarField Se = (SbSafe-Spc_irr_)/(Spc_max_-Spc_irr_);
    // volScalarField SeSafe = Foam::max(Se, SMALL);

    dpcds = ((pc0_ * c_ * (-scalar(1.0) + Swc_ + Sgr_)) / Foam::pow(Sb - Swc_ + VSMALL, 2))*
                            Foam::pow((scalar(1.0) - Sb - Sgr_)/(Sb - Swc_),c_-1); 

    // correct pc
    pc = pc0_ * Foam::pow(scalar(1.0) - Sb - Sgr_, c_) / (Sb - Swc_ + VSMALL); 
}

} // End namespace Foam
