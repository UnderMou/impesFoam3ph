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
    Foam::isothermModel::Henry

Description
    Henry adsorption model for three-phase flow.

    Reference:
    \verbatim
        Kalam, S., Abu-Khamsin, S. A., Kamal, M. S., & Patil, S. (2021).
        Surfactant adsorption isotherms: A review.
        ACS omega, 6(48), 32342-32348.
    \endverbatim

    The default model coefficients are
    \verbatim

        Henry
        {
            Kd    0.05;
        }

    \endverbatim

SourceFiles
    Henry.C

\*---------------------------------------------------------------------------*/

#ifndef Henry_H
#define Henry_H

#include "isothermModel.H"

namespace Foam
{

class Henry
:
    public isothermModel
{
    scalar Kd_;

public:
    TypeName("Henry");

    Henry(const dictionary& dict, foamAuxFields* aux);

    virtual void correct
    (
       
    ) const;
};

} // End namespace Foam

#endif
