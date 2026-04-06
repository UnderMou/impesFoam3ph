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
    Foam::isothermModel

Description
    This class supplies an isotherm model template, 
    and is the base class for a number of isotherm models implementation.

SourceFiles
    isothermModel.C

\*---------------------------------------------------------------------------*/

#ifndef isothermModel_H
#define isothermModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "foamAuxFields.H"

namespace Foam
{

class isothermModel
{
protected:
    const dictionary& dict_;
    foamAuxFields* aux_;

public:
    TypeName("isothermModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        isothermModel,
        dictionary,
        (
            const dictionary& dict,
            foamAuxFields* aux
        ),
        (dict, aux)
    );

    isothermModel(const dictionary& dict, foamAuxFields* aux)
    :
        dict_(dict),
        aux_(aux)
    {}

    virtual ~isothermModel() = default;

    // ---- Interface ----

    virtual void correct
    (

    ) const = 0;

    // Factory
    static autoPtr<isothermModel> New
    (
        const dictionary& transportProperties,
        foamAuxFields* aux
    );
};

} // End namespace Foam

#endif
