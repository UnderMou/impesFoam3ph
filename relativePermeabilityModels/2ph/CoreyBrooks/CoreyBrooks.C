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

#include "CoreyBrooks.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    CoreyBrooks,
    0
);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    CoreyBrooks,
    dictionary
);

CoreyBrooks::CoreyBrooks(const dictionary& dict)
:
    relativePermeabilityModel(dict),
    Sa_min_(readScalar(dict.lookup("Sa_min"))),
    Sb_min_(readScalar(dict.lookup("Sb_min"))),
    a_exp_(readScalar(dict.lookup("a_exp"))),
    b_exp_(readScalar(dict.lookup("b_exp"))),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_(readScalar(dict.lookup("krb_max")))
{}

void CoreyBrooks::correct
(
    volScalarField& kra,
    volScalarField& krb,
    volScalarField& Sb,
    volScalarField& Csw
) const
{

    // Effective saturations (normalization)
    volScalarField Se = (Sb - Sb_min_) / (1.0 - Sa_min_ - Sb_min_);

    // Check
    Se = min(max(Se, scalar(0.0)), scalar(1.0));

    // Corey-type relative permeabilities
    kra = kra_max_ * Foam::pow(scalar(1.0) - Se, a_exp_);
    krb = krb_max_ * Foam::pow(Se, b_exp_); 

}

} // End namespace Foam
