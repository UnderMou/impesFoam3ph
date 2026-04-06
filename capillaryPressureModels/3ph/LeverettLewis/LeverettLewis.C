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

#include "LeverettLewis.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    LeverettLewis,
    0
);

addToRunTimeSelectionTable
(
    capillaryPressureModel,
    LeverettLewis,
    dictionary
);

LeverettLewis::LeverettLewis(const dictionary& dict)
:
    capillaryPressureModel(dict),
    epsilon_(readScalar(dict.lookup("epsilon")))
{

}

void LeverettLewis::correct
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
    // correct pcwo and its derivatives
    pcow = -5.0 * epsilon_ * (scalar(2.0) - Sb) * (scalar(1.0) - Sb);
    dpcow_dsw = -5.0 * epsilon_ * (scalar(2.0)*Sb - scalar(3.0));
    dpcow_dsg = zeroField;

    // correct pcgo and its derivatives
    pcgo = 5.0 * epsilon_ * (scalar(2.0) - Sa) * (scalar(1.0) - Sa);
    dpcgo_dsg = 5.0 * epsilon_ * (scalar(2.0)*Sa - scalar(3.0));
    dpcgo_dsw = zeroField;
}

} // End namespace Foam
