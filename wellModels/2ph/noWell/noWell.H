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
    Foam::wellModel::noWell

Description
    no well model for three-phase flow fow when not considering injection wells.

    The default model coefficients are
    \verbatim

        wellModel noWell;

        noWell
        {
            
        }

    \endverbatim

SourceFiles
    noWell.C

\*---------------------------------------------------------------------------*/

#ifndef noWell_H
#define noWell_H

#include "wellModel.H"

namespace Foam
{

class noWell
:
    public wellModel
{

public:
    TypeName("noWell");

    noWell(const dictionary& dict);

    virtual void correct
    (
        volScalarField& qb,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod,
        scalar t
    ) const;
};

} // End namespace Foam

#endif
