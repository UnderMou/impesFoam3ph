/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "prescribedFracFlowVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prescribedFracFlowVelocity::
prescribedFracFlowVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(p.size(), vector::zero),
    fractionalFlow_(0)
{}


Foam::prescribedFracFlowVelocity::
prescribedFracFlowVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    U_("U", dict, p.size()),
    fractionalFlow_(dict.lookupOrDefault<scalar>("fractionalFlow", 0))
{}


// Foam::prescribedFracFlowVelocity::
// prescribedFracFlowVelocity
// (
//     const prescribedFracFlowVelocity& ptf,
//     const fvPatch& p,
//     const DimensionedField<vector, volMesh>& iF,
//     const fvPatchFieldMapper& mapper
// )
// :
//     fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
//     U_(ptf.U_, vector::zero),
//     fractionalFlow_(ptf.fractionalFlow_)
// {}


Foam::prescribedFracFlowVelocity::
prescribedFracFlowVelocity
(
    const prescribedFracFlowVelocity& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    U_(ptf.U_),
    fractionalFlow_(ptf.fractionalFlow_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::prescribedFracFlowVelocity::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(fractionalFlow_*U_);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::prescribedFracFlowVelocity::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    // os.writeKeyword("U")
    //     << U_ << token::END_STATEMENT << nl;
    os.writeKeyword("U")
        << "uniform " << U_[0]
        << token::END_STATEMENT << nl;

    os.writeKeyword("fractionalFlow")
        << fractionalFlow_ << token::END_STATEMENT << nl;
    
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       prescribedFracFlowVelocity
   );
}


// ************************************************************************* //
