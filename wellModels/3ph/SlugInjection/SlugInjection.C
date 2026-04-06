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

#include "SlugInjection.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    SlugInjection,
    0
);

addToRunTimeSelectionTable
(
    wellModel,
    SlugInjection,
    dictionary
);

SlugInjection::SlugInjection(const dictionary& dict)
:
    wellModel(dict),
    slugTime_w_(readScalar(dict.lookup("slugTime_w"))),
    Fb_inj_(readScalar(dict.lookup("Fb_inj")))
{
    scalar cycleTime_ = slugTime_w_ + slugTime_g_;
    scalar t_cycle_ = std::fmod(0, cycleTime_);
}

void SlugInjection::correct
(
    volScalarField& qa,
        volScalarField& qb,
        const volScalarField& Fa,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod,
        scalar t
) const
{

    // Map time into cycle
    t_cycle_ = std::fmod(t, cycleTime_);

    // Check which slug we are in
    bool inWater = (t_cycle_ >= 0.0 && t_cycle_ < slugTime_w_);
    bool inGas   = (t_cycle_ >= slugTime_w_  && t_cycle_ < cycleTime_);

    if inWater
    {
        qb = qt_inj - Fb * qt_prod; 
        qa = - Fa * qt_prod;
    }
    
    if inGas
    {
        qb = - Fb * qt_prod; 
        qa = qt_inj - Fa * qt_prod;
    }

}

} // End namespace Foam
