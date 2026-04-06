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
    Foam::surfactantTransportModel

Description
    This class supplies a surfactant transport equation model template, 
    and is the base class for a number of surfactant transport models implementation.

SourceFiles
    surfactantTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfactantTransportModel_H
#define surfactantTransportModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "foamAuxFields.H"
#include "isothermModel.H"

namespace Foam
{

class surfactantTransportModel
{
protected:
    const dictionary& dict_;
    foamAuxFields* aux_;
    autoPtr<isothermModel> isothermModel_;

    dimensionedScalar rho_b_;

public:
    TypeName("surfactantTransportModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        surfactantTransportModel,
        dictionary,
        (
            const dictionary& dict,
            foamAuxFields* aux
        ),
        (dict, aux)
    );

    surfactantTransportModel(const dictionary& dict, foamAuxFields* aux)
    :
        dict_(dict),
        aux_(aux),
        rho_b_(dict.parent().lookup("rho_b"))
    {
        isothermModel_ = isothermModel::New(dict, aux);
    }

    virtual ~surfactantTransportModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        const volScalarField& Sb,
        const surfaceScalarField& phib,
        const volScalarField& eps
    ) const = 0;

    // Factory
    static autoPtr<surfactantTransportModel> New
    (
        const dictionary& transportProperties,
        foamAuxFields* aux
    );
};

} // End namespace Foam

#endif
