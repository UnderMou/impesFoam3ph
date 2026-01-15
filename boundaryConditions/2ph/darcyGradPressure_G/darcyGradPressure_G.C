#include "darcyGradPressure_G.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyGradPressure_G::darcyGradPressure_G
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_("Mf"),
    phiName_("phi"),
    phiGfName_("phiG")
{}

Foam::darcyGradPressure_G::darcyGradPressure_G
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    phiGfName_(dict.lookupOrDefault<word>("phiG","phiG"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::darcyGradPressure_G::darcyGradPressure_G
(
    const darcyGradPressure_G& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_)
{}

// Foam::darcyGradPressure_G::darcyGradPressure_G
// (
//     const darcyGradPressure_G& ptf
// )
//     :
//     fixedGradientFvPatchScalarField(ptf),
//     MfName_(ptf.MfName_),
//     phiName_(ptf.phiName_)
// {}

Foam::darcyGradPressure_G::darcyGradPressure_G
(
    const darcyGradPressure_G& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyGradPressure_G::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& Mf=
        patch().lookupPatchField<surfaceScalarField, scalar>(MfName_);

    const fvsPatchField<scalar>& phi=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvsPatchField<scalar>& phiGf=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiGfName_);

    gradient() = - (phi-phiGf)/(Mf+ROOTVSMALL)/(patch().magSf());

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void Foam::darcyGradPressure_G::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "phiG", "phiG", phiGfName_);
    writeEntry(os, "Mf", MfName_);
    writeEntry(os, "phi", phiName_);
    writeEntry(os, "phiG", phiGfName_);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyGradPressure_G
);
}

// ************************************************************************* //
