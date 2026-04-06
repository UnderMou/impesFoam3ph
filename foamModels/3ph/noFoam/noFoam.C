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
    const volScalarField& eps,
    const volScalarField& K,
    const volVectorField& gradP
) const
{
    // no foam
    // pass
}

} // End namespace Foam
