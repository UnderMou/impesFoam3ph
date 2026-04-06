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
    Foam::wellModel

Description
    This class supplies a well model template, 
    and is the base class for a number of well models implementation.

SourceFiles
    wellModelModel.C

\*---------------------------------------------------------------------------*/

#ifndef wellModel_H
#define wellModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class wellModel
{
protected:
    const dictionary& dict_;

public:
    TypeName("wellModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        wellModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    wellModel(const dictionary& dict)
    :
        dict_(dict)
    {}

    virtual ~wellModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        volScalarField& qb,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod,
        scalar t
    ) const = 0;

    // Factory
    static autoPtr<wellModel> New
    (
        const dictionary& wellProperties
    );
};

} // End namespace Foam

#endif
