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

#include "vanGenuchten.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    vanGenuchten,
    0
);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    vanGenuchten,
    dictionary
);

vanGenuchten::vanGenuchten(const dictionary& dict)
:
    relativePermeabilityModel(dict),
    Sa_min_(readScalar(dict.lookup("Sa_min"))),
    Sb_min_(readScalar(dict.lookup("Sb_min"))),
    n_(readScalar(dict.lookup("n"))),
    m_(1.0 - 1.0/n_),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_(readScalar(dict.lookup("krb_max")))
{}

void vanGenuchten::correct
(
    volScalarField& kra,
    volScalarField& krb,
    volScalarField& Sb,
    volScalarField& Csw
) const
{

    // Calculate effective saturation (normalization)
    volScalarField Se = (Sb - Sb_min_) / (1.0 - Sa_min_ - Sb_min_);

    // Check
    Se = min(max(Se, scalar(0.0)), scalar(1.0));

    // Calculate relative permeabilities
    kra = kra_max_ * Foam::pow(scalar(1.0) - Se, 0.5) * Foam::pow(scalar(1.0) - Foam::pow(Se, 1.0/m_), 2*m_);
    krb = krb_max_ * Foam::pow(Se, 0.5) * Foam::pow(scalar(1.0) - Foam::pow(scalar(1.0) - Foam::pow(Se, 1.0/m_), m_), 2); 

}

} // End namespace Foam
