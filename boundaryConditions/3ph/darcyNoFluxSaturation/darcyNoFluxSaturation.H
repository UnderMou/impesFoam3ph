/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::darcyNoFluxSaturation

Description
    Coupled fixed-gradient boundary condition for the two independent
    saturations Sa and Sb of a three-phase Darcy-flow model.

    Compatible with the following solver flux definitions:

        phi   = phiP + phiG + phiPc_a - phiPc_b

        phia  = Faf*phiP + (Laf/Lf)*phiG + phiPc_a

        phib  = Fbf*phiP + (Lbf/Lf)*phiG - phiPc_b

    with

        phiPc_a =
            Maf*(dpcgo_dsw_f*snGrad(Sb)
               + dpcgo_dsg_f*snGrad(Sa))*magSf

        phiPc_b =
            Mbf*(dpcow_dsw_f*snGrad(Sb)
               + dpcow_dsg_f*snGrad(Sa))*magSf.

    Together with a pressure boundary condition imposing phi = 0,
    this boundary condition computes snGrad(Sa) and snGrad(Sb)
    so that phia = 0 and phib = 0.

    The same boundary-condition class is used for Sa and Sb.
    Set "component Sa" for Sa and "component Sb" for Sb.

    The condition automatically reduces from the full 2x2 three-phase
    system to a single independent phase no-flow equation when one phase
    becomes immobile. Phase activity is determined from Maf, Mbf and Mcf
    using a relative mobilityTolerance.

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
    // Private data

        word component_;

        word SaName_;
        word SbName_;

        word FafName_;
        word FbfName_;

        word LafName_;
        word LbfName_;
        word LfName_;

        word MafName_;
        word MbfName_;
        word McfName_;

        word phiGName_;

        word dpcgoDswName_;
        word dpcgoDsgName_;
        word dpcowDswName_;
        word dpcowDsgName_;

        scalar SaMin_;
        scalar SaMax_;
        scalar SbMin_;
        scalar SbMax_;
        scalar ScMin_;
        scalar ScMax_;

        bool limitSaturation_;

        //- Relative mobility threshold used to decide whether a phase is active
        scalar mobilityTolerance_;

        scalar capillaryTolerance_;
        scalar singularTolerance_;


public:

    TypeName("darcyNoFluxSaturation");


    // Constructors

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


    // Member functions

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
