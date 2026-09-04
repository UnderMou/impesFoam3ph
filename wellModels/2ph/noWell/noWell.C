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

List<well> noWell::readWells(const dictionary& wellsDict)
{
    return List<well>(0);
}

noWell::noWell(const dictionary& dict)
:
    wellModel(dict)
{
    wells_ = readWells(wellsDict_);
    // checkRateBalance();
}

void noWell::source_pEqn
(
    fvScalarMatrix& pEqn,
    const volScalarField& p,
    const volScalarField& mob_t,
    const volScalarField& WI,
    volScalarField& wellCoeff,
    volScalarField& wellSource,
    const scalar& rho_a,
    const scalar& rho_b,
    const volScalarField& mob_a,
    const volScalarField& mob_b,
    const dimensionedVector& g
) 
{
    // no well
}

void noWell::source_SbEqn
(
    fvScalarMatrix& SbEqn,
    const volScalarField& Sb, 
    const volScalarField& Fb, 
    const volScalarField& p,
    scalar t,
    volScalarField& qb
)
{
    // no well
}

void noWell::correct
(
    volScalarField& qt,
    volScalarField& qb,
    const volScalarField& Fb,
    const volScalarField& p, 
    scalar t,
    const volScalarField& mob_t,
    const volScalarField& WI,
    volScalarField& p_bh,
    volScalarField& qs,
    const volScalarField& Cs,
    const scalar& rho_a,
    const scalar& rho_b,
    const volScalarField& mob_a,
    const volScalarField& mob_b,
    const dimensionedVector& g
) 
{
    // no well
}

} // End namespace Foam
