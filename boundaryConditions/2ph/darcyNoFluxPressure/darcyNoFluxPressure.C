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

#include "darcyNoFluxPressure.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyNoFluxPressure::darcyNoFluxPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    MfName_("Mf"),
    phiGName_("phiG"),
    phiPcName_("phiPc")
{}


Foam::darcyNoFluxPressure::darcyNoFluxPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    // gradientRequired=false: the gradient is computed by this BC
    fixedGradientFvPatchScalarField(p, iF, dict, false),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    phiGName_(dict.lookupOrDefault<word>("phiG", "phiG")),
    phiPcName_(dict.lookupOrDefault<word>("phiPc", "phiPc"))
{
    // Preserve value/gradient when reading a written time directory.
    // Otherwise start from the adjacent-cell value with zero gradient;
    // updateCoeffs() will replace the gradient before matrix assembly.
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


Foam::darcyNoFluxPressure::darcyNoFluxPressure
(
    const darcyNoFluxPressure& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    MfName_(ptf.MfName_),
    phiGName_(ptf.phiGName_),
    phiPcName_(ptf.phiPcName_)
{}


Foam::darcyNoFluxPressure::darcyNoFluxPressure
(
    const darcyNoFluxPressure& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    phiGName_(ptf.phiGName_),
    phiPcName_(ptf.phiPcName_)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyNoFluxPressure::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& Mf =
        patch().lookupPatchField<surfaceScalarField, scalar>(MfName_);

    const fvsPatchField<scalar>& phiG =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiGName_);

    const fvsPatchField<scalar>& phiPc =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiPcName_);

    /*
        Total flux on the boundary:

            phi = phiP + phiG + phiPc

        Pressure contribution from

            fvm::laplacian(-Mf, p)

        is

            phiP = -Mf*snGrad(p)*magSf.

        Impermeability requires

            phi = 0

        therefore

            -Mf*snGrad(p)*magSf + phiG + phiPc = 0

        and hence

            snGrad(p) = (phiG + phiPc)/(Mf*magSf).
    */

    gradient() = (phiG + phiPc)/((Mf + ROOTVSMALL)*patch().magSf());

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::darcyNoFluxPressure::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "phiG", "phiG", phiGName_);
    writeEntryIfDifferent<word>(os, "phiPc", "phiPc", phiPcName_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyNoFluxPressure
);
}

// ************************************************************************* //
