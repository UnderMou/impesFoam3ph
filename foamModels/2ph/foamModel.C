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

#include "foamModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    foamModel,
    0
);

defineRunTimeSelectionTable
(
    foamModel,
    dictionary
);

autoPtr<foamModel>
foamModel::New
(
    const dictionary& transportProperties,
    foamAuxFields* aux
)
{
    word modelName
    (
        transportProperties.lookup("foamModel")
    );

    Info<< "Selecting foam model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        transportProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown foam model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<foamModel>
    (
        iter()(modelDict, aux)
    );
}


} // End namespace Foam
