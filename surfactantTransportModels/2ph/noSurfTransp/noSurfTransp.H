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
    Foam::surfactantTransportModel::noSurfTransp

Description
    No tansport of surfactant model for when it is neglected in two-phase flow.

    The default model coefficients are
    \verbatim

        noSurfTransp
        {

        }

    \endverbatim

SourceFiles
    noSurfTransp.C

\*---------------------------------------------------------------------------*/

#ifndef noSurfTransp_H
#define noSurfTransp_H

#include "fvCFD.H"
#include "surfactantTransportModel.H"

namespace Foam
{

class noSurfTransp
:
    public surfactantTransportModel
{

public:
    TypeName("noSurfTransp");

    noSurfTransp(const dictionary& dict, foamAuxFields* aux);

    // noSurfTransp functions

    virtual void correct
    (
        const volScalarField& Sb,
        const surfaceScalarField& phib,
        const volScalarField& eps
    ) const;
};

} // End namespace Foam

#endif
