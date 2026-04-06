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
    Foam::isothermModel::noAdsorption

Description
    no adsorption model for three-phase flow for when adsorption effects 
    are neglected.


    The default model coefficients are
    \verbatim

        noAdsorption
        {

        }

    \endverbatim

SourceFiles
    noAdsorption.C

\*---------------------------------------------------------------------------*/

#ifndef noAdsorption_H
#define noAdsorption_H

#include "isothermModel.H"

namespace Foam
{

class noAdsorption
:
    public isothermModel
{

public:
    TypeName("noAdsorption");

    noAdsorption(const dictionary& dict, foamAuxFields* aux);

    virtual void correct
    (
       
    ) const;
};

} // End namespace Foam

#endif
