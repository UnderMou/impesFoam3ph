#include "Csw.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    Csw,
    0
);

addToRunTimeSelectionTable
(
    surfactantTransportModel,
    Csw,
    dictionary
);

Csw::Csw(const dictionary& dict, foamAuxFields* aux)
:
    surfactantTransportModel(dict, aux),
    rho_s_(readScalar(dict.lookup("rho_s")))
{
    rho_sw_ = rho_s_ / rho_b_.value();
    Info << "rho_sw = " << rho_sw_ << endl;
}

void Csw::correct
(
    const volScalarField& Sb,
    const surfaceScalarField& phib,
    const volScalarField& eps
) const
{
    if (!aux_ || !aux_->Cs)
    {
        FatalErrorInFunction
            << "Surfactant transport model requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& Cs        = *aux_->Cs;
    volScalarField& dCsEqdCs  = *aux_->dCsEqdCs;
    volScalarField& Fcsw      = *aux_->Fcsw;
    volScalarField& AcumCoeff = *aux_->AcumCoeff_csw;

    Fcsw = -fvc::div(phib); // + qb
    AcumCoeff = eps*Sb + (scalar(1.0) - eps)*rho_sw_*dCsEqdCs;

    fvScalarMatrix CsEqn
    (
        AcumCoeff*fvm::ddt(Cs) + fvc::div(phib,Cs) 
        ==
        fvc::Sp(-Fcsw,Cs)
    );

    CsEqn.solve();

    Info << "Surfactant concentration: " << " Min(Cs) = " << gMin(Cs) << " Max(Cs) = " << gMax(Cs) << endl;
}

} // End namespace Foam
