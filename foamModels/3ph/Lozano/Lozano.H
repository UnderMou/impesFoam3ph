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
    Foam::foamModel::Lozano

Description
    Lozano foam model for three-phase flow.

    Reference:
    \verbatim
        Lozano, L. F., Chapiro, G., & Marchesin, D. (2025). 
        The Riemann problem for three-phase foam flow in porous media. 
        arXiv preprint arXiv:2506.08152.
    \endverbatim

    The default model coefficients are
    \verbatim

        Lozano
        {
            MRF     1750;
        }

    \endverbatim

SourceFiles
    Lozano.C

\*---------------------------------------------------------------------------*/

#ifndef Lozano_H
#define Lozano_H

#include "foamModel.H"

namespace Foam
{

class Lozano
:
    public foamModel
{
    scalar MRF_;

public:
    TypeName("Lozano");

    Lozano(const dictionary& dict, foamAuxFields* aux);

    // Correct
    virtual void correct
    (
        volScalarField& kra,
        const volVectorField U,
        const volScalarField& Sa,
        const volScalarField& Sb,
        volScalarField& Sc,
        const surfaceScalarField& phia,
        const volScalarField& eps,
        const volScalarField& K,
        const volVectorField& gradP
    ) const;
};

} // End namespace Foam

#endif
