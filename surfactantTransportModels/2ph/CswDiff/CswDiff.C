/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CswDiff.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    CswDiff,
    0
);

addToRunTimeSelectionTable
(
    surfactantTransportModel,
    CswDiff,
    dictionary
);

CswDiff::CswDiff(const dictionary& dict, foamAuxFields* aux)
:
    surfactantTransportModel(dict, aux),
    rho_s_(readScalar(dict.lookup("rho_s"))),
    D_(readScalar(dict.lookup("D"))),
    Deff_("Deff", dimArea/dimTime, D_)
{
    rho_sw_ = rho_s_ / rho_b_.value();
    Info << "rho_sw = " << rho_sw_ << endl;
    Info << "Deff = " << Deff_ << endl;
}

void CswDiff::correct
(
    const volScalarField& Sb,
    const surfaceScalarField& phib,
    const volScalarField& eps,
    const volScalarField& qb,
    const volScalarField& qs
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

    // isotherm
    Info<< "Using isotherm model: " << isothermModel_->type() << nl << endl;
    isothermModel_->correct();

    Fcsw = qb - fvc::div(phib);
    AcumCoeff = eps*Sb + (scalar(1.0) - eps)*rho_sw_*dCsEqdCs;

    // solve Cs equation
    fvScalarMatrix CsEqn
    (
        AcumCoeff*fvm::ddt(Cs) + fvc::div(phib,Cs) - fvm::laplacian(Deff_,Cs)
        ==
        // -Fcsw*Cs + qs
        fvc::Sp(-Fcsw,Cs) + qs
        // fvm::Sp(-Fcsw,Cs) + qs
    );
    CsEqn.solve();

    Info << "Surfactant concentration: " << " Min(Cs) = " << gMin(Cs) << " Max(Cs) = " << gMax(Cs) << endl;
}

} // End namespace Foam
