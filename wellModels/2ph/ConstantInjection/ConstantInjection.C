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

#include "ConstantInjection.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    ConstantInjection,
    0
);

addToRunTimeSelectionTable
(
    wellModel,
    ConstantInjection,
    dictionary
);

ConstantInjection::ConstantInjection(const dictionary& dict)
:
    wellModel(dict),
    Fb_inj_(readScalar(dict.lookup("Fb_inj")))
{

}

tmp<volScalarField> ConstantInjection::explicitSource_pEqn
(
    volScalarField& qb,
    const volScalarField& Sb, 
    const volScalarField& Fb, 
    const volScalarField& p,
    const volScalarField& qt_inj,
    const volScalarField& qt_prod,
    scalar t,
    const volScalarField& mob_t,
    const volScalarField& WI,
    const volScalarField& p_bh,
    const volScalarField& isWell
) const
{
    return qt_inj - qt_prod;
}

tmp<volScalarField> ConstantInjection::explicitSource_SEqn
(
    volScalarField& qb,
    const volScalarField& Sb, 
    const volScalarField& Fb, 
    const volScalarField& p,
    const volScalarField& qt_inj,
    const volScalarField& qt_prod,
    scalar t
) const
{
    // Calculate injection/production rates based on:
    //  > specified injection rates;
    //  > production rates as function of injection rates and fractional flow
    qb = Fb_inj_ * qt_inj - Fb * qt_prod; 
    return qb;
}

} // End namespace Foam
