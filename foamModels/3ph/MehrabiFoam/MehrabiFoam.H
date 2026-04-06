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
    Foam::foamModel::MehrabiFoam

Description
    Mehrabi foam model for three-phase flow.
    Similar to Foam::foamModel::Lozano, but it considers MRF = MRF(Cs)

    Reference:
    \verbatim
        Mehrabi, M., Sepehrnoori, K., & Delshad, M. (2022).
        Displacement theory of low-tension gas flooding.
        Transport in Porous Media, 142(3), 475-491.
    \endverbatim

    The default model coefficients are
    \verbatim

        MehrabiFoam
        {
        Cr       1e-2;
        MRF      10;
        }

    \endverbatim

SourceFiles
    MehrabiFoam.C

\*---------------------------------------------------------------------------*/

#ifndef MehrabiFoam_H
#define MehrabiFoam_H

#include "foamModel.H"

namespace Foam
{

class MehrabiFoam
:
    public foamModel
{
    scalar MRF_;
    scalar Cr_;

public:
    TypeName("MehrabiFoam");

    MehrabiFoam(const dictionary& dict, foamAuxFields* aux);

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
