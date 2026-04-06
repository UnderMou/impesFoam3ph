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
    Foam::foamModel

Description
    This class supplies a foam model template, 
    and is the base class for a number of foam models implementation.

SourceFiles
    foamModel.C

\*---------------------------------------------------------------------------*/

#ifndef foamModel_H
#define foamModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "foamAuxFields.H"

namespace Foam
{

class foamModel
{
protected:
    const dictionary& dict_;
    foamAuxFields* aux_;

    dimensionedScalar mu_b_;

public:
    TypeName("foamModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        foamModel,
        dictionary,
        (
            const dictionary& dict,
            foamAuxFields* aux
        ),
        (dict, aux)
    );

    foamModel(const dictionary& dict, foamAuxFields* aux)
    :
        dict_(dict),
        aux_(aux),
        mu_b_(dict.parent().lookup("mu_b"))
    {}

    virtual ~foamModel() = default;

    // ---- Interface ----
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
    ) const = 0;

    // Factory
    static autoPtr<foamModel> New
    (
        const dictionary& transportProperties,
        foamAuxFields* aux
    );
};

} // End namespace Foam

#endif
