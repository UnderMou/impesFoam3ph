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
    Foam::darcyNoFluxSaturation

Description
    Saturation boundary condition for an impermeable wall in two-phase Darcy
    flow.

    The phase-b face flux is assumed to be

        phib = Fbf*phiP + (Lbf/Lf)*phiG + phiPc

    and the pressure boundary condition imposes

        phi = phiP + phiG + phiPc = 0.

    With capillary effects,

        phiPc = Mbf*dpcdsf*snGrad(Sb)*magSf,

    and imposing phib = 0 gives

        snGrad(Sb)
          = (Fbf - Lbf/Lf)*phiG
           /[(1 - Fbf)*Mbf*dpcdsf*magSf].

    In the absence of a usable capillary coefficient, the saturation gradient
    cannot control phib. In that case this patch field uses a neutral zero
    gradient and the solver must impose phib = 0 explicitly on this patch.

Usage

    walls
    {
        type            darcyNoFluxSaturation;

        Fbf             Fbf;
        Lbf             Lbf;
        Lf              Lf;
        phiG            phiG;
        Mbf             Mbf;
        dpcdsf          dpcdsf;

        capillaryTolerance  1e-30;

        value           uniform 0.2;
    }

SourceFiles
    darcyNoFluxSaturation.C

\*---------------------------------------------------------------------------*/

#ifndef darcyNoFluxSaturation_H
#define darcyNoFluxSaturation_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class darcyNoFluxSaturation
:
    public fixedGradientFvPatchScalarField
{
    word FbfName_;
    word LbfName_;
    word LfName_;
    word phiGName_;
    word MbfName_;
    word dpcdsfName_;

    scalar capillaryTolerance_;

public:

    TypeName("darcyNoFluxSaturation");

    darcyNoFluxSaturation
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    darcyNoFluxSaturation
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    darcyNoFluxSaturation
    (
        const darcyNoFluxSaturation&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    darcyNoFluxSaturation
    (
        const darcyNoFluxSaturation&
    ) = delete;

    darcyNoFluxSaturation
    (
        const darcyNoFluxSaturation&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new darcyNoFluxSaturation(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
