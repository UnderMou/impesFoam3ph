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

Class
    Foam::wellModel::SlugInjection

Description
    Constant injection rate model (water/gas) given a constant time interval for three-phase flow.

    The default model coefficients are
    \verbatim

        wellModel SlugInjection;

        SlugInjection
        {
            slugTime_w;
            slugTime_g;            
        }

    \endverbatim

SourceFiles
    SlugInjection.C

\*---------------------------------------------------------------------------*/

#ifndef SlugInjection_H
#define SlugInjection_H

#include "wellModel.H"

namespace Foam
{

class SlugInjection
:
    public wellModel
{
    scalar slugTime_w_;
    scalar slugTime_g_;
    scalar cycleTime_;
    scalar t_cycle_;
    bool inWater_; 
    bool inGas_; 

public:
    TypeName("SlugInjection");

    SlugInjection(const dictionary& dict);

    virtual void correct
    (
        volScalarField& qa,
        volScalarField& qb,
        const volScalarField& Fa,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod,
        scalar t
    ) const;
};

} // End namespace Foam

#endif
