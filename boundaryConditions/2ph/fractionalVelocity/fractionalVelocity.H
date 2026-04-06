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
    Foam::fractionalVelocity

Description
    This boundary condition provides a velocity boundary condition for
    multiphase flow based on calculated phase fractional flow and
    a total darcy velocity.

Usage
    \table
        Property        | Description                  | Required    | Default value
        U               | total darcy velocity [m/s]   | yes         |
        fractionalFlow  | phase fracitonal flow        | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fractionalVelocity;
        U               (1e-4 0 0);
        fractionalFlow  fi; // i-th phase fractional flow; i = a,b
        value           uniform (0 0 0); // placeholder
    }
    \endverbatim

SourceFiles
    fractionalVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef fractionalVelocity_H
#define fractionalVelocity_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
    Class fractionalVelocity Declaration
\*---------------------------------------------------------------------------*/

class fractionalVelocity
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Total darcy velocity
        vectorField U_;

        //- Name of the phase fractional flow
        word fractionalFlow_;


public:

   //- Runtime type information
   TypeName("fractionalVelocity");


   // Constructors

        //- Construct from patch and internal field
        fractionalVelocity
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fractionalVelocity
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fractionalVelocity
        //  onto a new patch
        fractionalVelocity
        (
            const fractionalVelocity&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fractionalVelocity
        (
            const fractionalVelocity&
        ) = delete;

        //- Copy constructor setting internal field reference
        fractionalVelocity
        (
            const fractionalVelocity&,
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
                new fractionalVelocity
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
