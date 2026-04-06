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

\*---------------------------------------------------------------------------*/

#include "relativePermeabilityModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    relativePermeabilityModel,
    0
);

defineRunTimeSelectionTable
(
    relativePermeabilityModel,
    dictionary
);

autoPtr<relativePermeabilityModel>
relativePermeabilityModel::New
(
    const dictionary& transportProperties
)
{
    word modelName
    (
        transportProperties.lookup("relativePermeabilityModel")
    );

    Info<< "Selecting relative permeability model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        transportProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown relative permeability model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<relativePermeabilityModel>
    (
        iter()(modelDict)
    );
}

} // End namespace Foam
