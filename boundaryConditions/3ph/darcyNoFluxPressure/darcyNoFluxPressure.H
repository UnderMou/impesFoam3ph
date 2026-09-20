/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Class
    Foam::darcyNoFluxPressure

Description
    Fixed-gradient pressure boundary condition for an impermeable Darcy-flow
    boundary when the total face flux is

        phi = phiP + phiG + phiPc

    with

        phiP = -Mf*snGrad(p)*magSf.

    The condition imposes phi = 0 by setting

        snGrad(p) = (phiG + phiPc)/(Mf*magSf).

    Thus gravity and capillarity are not forced individually to zero; the
    pressure-gradient contribution balances them.

Usage

    <patchName>
    {
        type            darcyNoFluxPressure;

        // Optional field names:
        Mf              Mf;
        phiG            phiG;
        phiPc           phiPc;

        value           uniform 0;
    }

SourceFiles
    darcyNoFluxPressure.C

\*---------------------------------------------------------------------------*/

#ifndef darcyNoFluxPressure_H
#define darcyNoFluxPressure_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class darcyNoFluxPressure
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of total face mobility field
        word MfName_;

        //- Name of gravitational face-flux field
        word phiGName_;

        //- Name of capillary face-flux field
        word phiPcName_;


public:

    //- Runtime type information
    TypeName("darcyNoFluxPressure");


    // Constructors

        //- Construct from patch and internal field
        darcyNoFluxPressure
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        darcyNoFluxPressure
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        darcyNoFluxPressure
        (
            const darcyNoFluxPressure&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without resetting the internal-field reference
        darcyNoFluxPressure
        (
            const darcyNoFluxPressure&
        ) = delete;

        //- Copy constructor setting internal-field reference
        darcyNoFluxPressure
        (
            const darcyNoFluxPressure&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal-field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new darcyNoFluxPressure(*this, iF)
            );
        }


    // Member functions

        //- Update the pressure gradient so that total flux is zero
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
