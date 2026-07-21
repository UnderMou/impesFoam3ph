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

#include "fractionalVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fractionalVelocity::
fractionalVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(p.size(), vector::zero),
    fractionalFlow_("none")
{}


Foam::fractionalVelocity::
fractionalVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    U_("U", dict, p.size()),
    fractionalFlow_(dict.lookupOrDefault<word>("fractionalFlow", "none"))
{}


// Foam::fractionalVelocity::
// fractionalVelocity
// (
//     const fractionalVelocity& ptf,
//     const fvPatch& p,
//     const DimensionedField<vector, volMesh>& iF,
//     const fvPatchFieldMapper& mapper
// )
// :
//     fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
//     U_(ptf.U_, vector::zero),
//     fractionalFlow_(ptf.fractionalFlow_)
// {}


Foam::fractionalVelocity::
fractionalVelocity
(
    const fractionalVelocity& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    U_(ptf.U_),
    fractionalFlow_(ptf.fractionalFlow_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fractionalVelocity::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const word fluxName = fractionalFlow_ + "f";

    const fvsPatchScalarField& fS =
        patch().lookupPatchField<surfaceScalarField, scalar>(fluxName);

    operator==(fS*U_);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::fractionalVelocity::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    // os.writeKeyword("U")
        // << U_ << token::END_STATEMENT << nl;
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
       fractionalVelocity
   );
}


// ************************************************************************* //
