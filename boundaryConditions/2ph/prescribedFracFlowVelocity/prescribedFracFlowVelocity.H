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
    Foam::prescribedFracFlowVelocity

Description
    This boundary condition provides a velocity boundary condition for
    multiphase flow based on prescribed phase fractional flow and
    a total darcy velocity.

Usage
    \table
        Property        | Description                             | Required    | Default value
        U               | total darcy velocity [m/s]              | yes         |
        fractionalFlow  | prescribed phase fracitonal flow        | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            prescribedFracFlowVelocity;
        U               uniform (1e-4 0 0);
        fractionalFlow  0.7; // phase fractional flow
        value           uniform (0 0 0); // placeholder
    }
    \endverbatim

SourceFiles
    prescribedFracFlowVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedFracFlowVelocity_H
#define prescribedFracFlowVelocity_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
    Class prescribedFracFlowVelocity Declaration
\*---------------------------------------------------------------------------*/

class prescribedFracFlowVelocity
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Total darcy velocity
        vectorField U_;

        //- Prescribed phase fractional flow
        scalar fractionalFlow_;


public:

   //- Runtime type information
   TypeName("prescribedFracFlowVelocity");


   // Constructors

        //- Construct from patch and internal field
        prescribedFracFlowVelocity
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        prescribedFracFlowVelocity
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  prescribedFracFlowVelocity
        //  onto a new patch
        prescribedFracFlowVelocity
        (
            const prescribedFracFlowVelocity&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        prescribedFracFlowVelocity
        (
            const prescribedFracFlowVelocity&
        ) = delete;

        //- Copy constructor setting internal field reference
        prescribedFracFlowVelocity
        (
            const prescribedFracFlowVelocity&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new prescribedFracFlowVelocity
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
