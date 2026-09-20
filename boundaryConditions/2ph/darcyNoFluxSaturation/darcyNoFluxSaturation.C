/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "darcyNoFluxSaturation.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

Foam::darcyNoFluxSaturation::darcyNoFluxSaturation
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    FbfName_("Fbf"),
    LbfName_("Lbf"),
    LfName_("Lf"),
    phiGName_("phiG"),
    MbfName_("Mbf"),
    dpcdsfName_("dpcdsf"),
    capillaryTolerance_(SMALL)
{}


Foam::darcyNoFluxSaturation::darcyNoFluxSaturation
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict, false),
    FbfName_(dict.lookupOrDefault<word>("Fbf", "Fbf")),
    LbfName_(dict.lookupOrDefault<word>("Lbf", "Lbf")),
    LfName_(dict.lookupOrDefault<word>("Lf", "Lf")),
    phiGName_(dict.lookupOrDefault<word>("phiG", "phiG")),
    MbfName_(dict.lookupOrDefault<word>("Mbf", "Mbf")),
    dpcdsfName_(dict.lookupOrDefault<word>("dpcdsf", "dpcdsf")),
    capillaryTolerance_
    (
        dict.lookupOrDefault<scalar>("capillaryTolerance", SMALL)
    )
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );

        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = Zero;
    }
}


Foam::darcyNoFluxSaturation::darcyNoFluxSaturation
(
    const darcyNoFluxSaturation& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    FbfName_(ptf.FbfName_),
    LbfName_(ptf.LbfName_),
    LfName_(ptf.LfName_),
    phiGName_(ptf.phiGName_),
    MbfName_(ptf.MbfName_),
    dpcdsfName_(ptf.dpcdsfName_),
    capillaryTolerance_(ptf.capillaryTolerance_)
{}


Foam::darcyNoFluxSaturation::darcyNoFluxSaturation
(
    const darcyNoFluxSaturation& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    FbfName_(ptf.FbfName_),
    LbfName_(ptf.LbfName_),
    LfName_(ptf.LfName_),
    phiGName_(ptf.phiGName_),
    MbfName_(ptf.MbfName_),
    dpcdsfName_(ptf.dpcdsfName_),
    capillaryTolerance_(ptf.capillaryTolerance_)
{}


void Foam::darcyNoFluxSaturation::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& Fbf =
        patch().lookupPatchField<surfaceScalarField, scalar>(FbfName_);

    const fvsPatchField<scalar>& Lbf =
        patch().lookupPatchField<surfaceScalarField, scalar>(LbfName_);

    const fvsPatchField<scalar>& Lf =
        patch().lookupPatchField<surfaceScalarField, scalar>(LfName_);

    const fvsPatchField<scalar>& phiG =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiGName_);

    const fvsPatchField<scalar>& Mbf =
        patch().lookupPatchField<surfaceScalarField, scalar>(MbfName_);

    const fvsPatchField<scalar>& dpcdsf =
        patch().lookupPatchField<surfaceScalarField, scalar>(dpcdsfName_);

    scalarField& gradSb = gradient();

    forAll(gradSb, facei)
    {
        const scalar Fb = Fbf[facei];
        const scalar Lt = Lf[facei];
        const scalar Mb = Mbf[facei];
        const scalar dpc = dpcdsf[facei];

        // If capillary control is unavailable, snGrad(Sb) cannot control
        // the phase flux. Use a neutral gradient here; the solver must set
        // phib = 0 directly on this patch.
        if
        (
            mag(Mb*dpc) <= capillaryTolerance_
         || mag(1.0 - Fb) <= SMALL
         || mag(Lt) <= SMALL
        )
        {
            gradSb[facei] = 0.0;
        }
        else
        {
            const scalar betaB = Lbf[facei]/Lt;

            // From:
            //   phib = Fbf*phiP + betaB*phiG + phiPc
            //   phiP = -(phiG + phiPc)
            // gives:
            //   phib = (betaB - Fbf)*phiG + (1 - Fbf)*phiPc
            //
            // Setting phib = 0:
            //   phiPc = (Fbf - betaB)*phiG/(1 - Fbf)
            //
            // With:
            //   phiPc = Mbf*dpcdsf*snGrad(Sb)*magSf

            gradSb[facei] =
                (Fb - betaB)*phiG[facei]
               /
                (
                    (1.0 - Fb)
                   *Mb
                   *dpc
                   *patch().magSf()[facei]
                );
        }
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::darcyNoFluxSaturation::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "Fbf", "Fbf", FbfName_);
    writeEntryIfDifferent<word>(os, "Lbf", "Lbf", LbfName_);
    writeEntryIfDifferent<word>(os, "Lf", "Lf", LfName_);
    writeEntryIfDifferent<word>(os, "phiG", "phiG", phiGName_);
    writeEntryIfDifferent<word>(os, "Mbf", "Mbf", MbfName_);
    writeEntryIfDifferent<word>(os, "dpcdsf", "dpcdsf", dpcdsfName_);

    os.writeKeyword("capillaryTolerance")
        << capillaryTolerance_ << token::END_STATEMENT << nl;

    writeEntry(os, "value", *this);
}


namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyNoFluxSaturation
);
}

// ************************************************************************* //
