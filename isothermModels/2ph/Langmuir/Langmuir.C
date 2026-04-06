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

#include "Langmuir.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    Langmuir,
    0
);

addToRunTimeSelectionTable
(
    isothermModel,
    Langmuir,
    dictionary
);

Langmuir::Langmuir(const dictionary& dict, foamAuxFields* aux)
:
    isothermModel(dict, aux),
    qo_(readScalar(dict.lookup("qo"))),
    KL_(readScalar(dict.lookup("KL")))
{}

void Langmuir::correct
(
    
) const
{

    if (!aux_ || !aux_->dCsEqdCs)
    {
        FatalErrorInFunction
            << "Isother model requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& Cs              = *aux_->Cs;
    volScalarField& dCsEqdCs        = *aux_->dCsEqdCs;
    volScalarField& CsEq            = *aux_->CsEq;

    // update adsorbed concentration and its derivative
    dCsEqdCs = qo_ * KL_ / (Foam::pow(scalar(1.0) + KL_ * Cs, 2));
    CsEq = qo_ * KL_ * Cs / (scalar(1.0) + KL_ * Cs);

}

} // End namespace Foam
