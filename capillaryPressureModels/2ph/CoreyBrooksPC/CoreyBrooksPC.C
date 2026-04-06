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
    // correct pc and its derivatives
    volScalarField SbSafe = Foam::max(Sb, Spc_irr_ + SMALL);
    volScalarField Se = (SbSafe-Spc_irr_)/(Spc_max_-Spc_irr_);
    volScalarField SeSafe = Foam::max(Se, SMALL);
    dpcds = -alpha_ * pc0_ * Foam::pow(SeSafe, -alpha_ - 1.0) / (Spc_max_ - Spc_irr_);
    pc = pc0_ * Foam::pow(SeSafe, -alpha_);
}

} // End namespace Foam
