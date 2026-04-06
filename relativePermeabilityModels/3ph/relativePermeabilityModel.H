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
    Foam::relativePermeabilityModel

Description
    This class supplies a relative permeability model template, 
    and is the base class for a number of relative permeabilityModel.

SourceFiles
    relativePermeabilityModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativePermeabilityModel_H
#define relativePermeabilityModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class relativePermeabilityModel
{
protected:
    const dictionary& dict_;

public:
    TypeName("relativePermeabilityModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        relativePermeabilityModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    relativePermeabilityModel(const dictionary& dict)
    :
        dict_(dict)
    {}

    virtual ~relativePermeabilityModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        volScalarField& kra,
        volScalarField& krb,
        volScalarField& krc,
        volScalarField& Sa,
        volScalarField& Sb,
        volScalarField& Csw
    ) const = 0;

    // Factory
    static autoPtr<relativePermeabilityModel> New
    (
        const dictionary& transportProperties
    );
};

} // End namespace Foam

#endif
