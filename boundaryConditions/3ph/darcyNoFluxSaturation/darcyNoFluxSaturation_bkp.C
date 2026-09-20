/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "darcyNoFluxSaturation.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyNoFluxSaturation::darcyNoFluxSaturation
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    component_(iF.name()),
    SaName_("Sa"),
    SbName_("Sb"),
    FafName_("Faf"),
    FbfName_("Fbf"),
    LafName_("Laf"),
    LbfName_("Lbf"),
    LfName_("Lf"),
    MafName_("Maf"),
    MbfName_("Mbf"),
    phiGName_("phiG"),
    dpcgoDswName_("dpcgo_dsw_f"),
    dpcgoDsgName_("dpcgo_dsg_f"),
    dpcowDswName_("dpcow_dsw_f"),
    dpcowDsgName_("dpcow_dsg_f"),
    SaMin_(0.0),
    SaMax_(1.0),
    SbMin_(0.0),
    SbMax_(1.0),
    ScMin_(0.0),
    ScMax_(1.0),
    limitSaturation_(true),
    capillaryTolerance_(SMALL),
    singularTolerance_(1e-12)
{}


Foam::darcyNoFluxSaturation::darcyNoFluxSaturation
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict, false),
    component_(dict.lookupOrDefault<word>("component", iF.name())),
    SaName_(dict.lookupOrDefault<word>("Sa", "Sa")),
    SbName_(dict.lookupOrDefault<word>("Sb", "Sb")),
    FafName_(dict.lookupOrDefault<word>("Faf", "Faf")),
    FbfName_(dict.lookupOrDefault<word>("Fbf", "Fbf")),
    LafName_(dict.lookupOrDefault<word>("Laf", "Laf")),
    LbfName_(dict.lookupOrDefault<word>("Lbf", "Lbf")),
    LfName_(dict.lookupOrDefault<word>("Lf", "Lf")),
    MafName_(dict.lookupOrDefault<word>("Maf", "Maf")),
    MbfName_(dict.lookupOrDefault<word>("Mbf", "Mbf")),
    phiGName_(dict.lookupOrDefault<word>("phiG", "phiG")),
    dpcgoDswName_
    (
        dict.lookupOrDefault<word>("dpcgo_dsw", "dpcgo_dsw_f")
    ),
    dpcgoDsgName_
    (
        dict.lookupOrDefault<word>("dpcgo_dsg", "dpcgo_dsg_f")
    ),
    dpcowDswName_
    (
        dict.lookupOrDefault<word>("dpcow_dsw", "dpcow_dsw_f")
    ),
    dpcowDsgName_
    (
        dict.lookupOrDefault<word>("dpcow_dsg", "dpcow_dsg_f")
    ),
    SaMin_(dict.lookupOrDefault<scalar>("SaMin", 0.0)),
    SaMax_(dict.lookupOrDefault<scalar>("SaMax", 1.0)),
    SbMin_(dict.lookupOrDefault<scalar>("SbMin", 0.0)),
    SbMax_(dict.lookupOrDefault<scalar>("SbMax", 1.0)),
    ScMin_(dict.lookupOrDefault<scalar>("ScMin", 0.0)),
    ScMax_(dict.lookupOrDefault<scalar>("ScMax", 1.0)),
    limitSaturation_(dict.lookupOrDefault<bool>("limitSaturation", true)),
    capillaryTolerance_
    (
        dict.lookupOrDefault<scalar>("capillaryTolerance", SMALL)
    ),
    singularTolerance_
    (
        dict.lookupOrDefault<scalar>("singularTolerance", 1e-12)
    )
{
    if (component_ != "Sa" && component_ != "Sb")
    {
        FatalIOErrorInFunction(dict)
            << "component must be Sa or Sb, but found "
            << component_ << exit(FatalIOError);
    }

    if
    (
        SaMin_ > SaMax_
     || SbMin_ > SbMax_
     || ScMin_ > ScMax_
    )
    {
        FatalIOErrorInFunction(dict)
            << "Invalid saturation bounds."
            << exit(FatalIOError);
    }

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
    component_(ptf.component_),
    SaName_(ptf.SaName_),
    SbName_(ptf.SbName_),
    FafName_(ptf.FafName_),
    FbfName_(ptf.FbfName_),
    LafName_(ptf.LafName_),
    LbfName_(ptf.LbfName_),
    LfName_(ptf.LfName_),
    MafName_(ptf.MafName_),
    MbfName_(ptf.MbfName_),
    phiGName_(ptf.phiGName_),
    dpcgoDswName_(ptf.dpcgoDswName_),
    dpcgoDsgName_(ptf.dpcgoDsgName_),
    dpcowDswName_(ptf.dpcowDswName_),
    dpcowDsgName_(ptf.dpcowDsgName_),
    SaMin_(ptf.SaMin_),
    SaMax_(ptf.SaMax_),
    SbMin_(ptf.SbMin_),
    SbMax_(ptf.SbMax_),
    ScMin_(ptf.ScMin_),
    ScMax_(ptf.ScMax_),
    limitSaturation_(ptf.limitSaturation_),
    capillaryTolerance_(ptf.capillaryTolerance_),
    singularTolerance_(ptf.singularTolerance_)
{}


Foam::darcyNoFluxSaturation::darcyNoFluxSaturation
(
    const darcyNoFluxSaturation& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    component_(ptf.component_),
    SaName_(ptf.SaName_),
    SbName_(ptf.SbName_),
    FafName_(ptf.FafName_),
    FbfName_(ptf.FbfName_),
    LafName_(ptf.LafName_),
    LbfName_(ptf.LbfName_),
    LfName_(ptf.LfName_),
    MafName_(ptf.MafName_),
    MbfName_(ptf.MbfName_),
    phiGName_(ptf.phiGName_),
    dpcgoDswName_(ptf.dpcgoDswName_),
    dpcgoDsgName_(ptf.dpcgoDsgName_),
    dpcowDswName_(ptf.dpcowDswName_),
    dpcowDsgName_(ptf.dpcowDsgName_),
    SaMin_(ptf.SaMin_),
    SaMax_(ptf.SaMax_),
    SbMin_(ptf.SbMin_),
    SbMax_(ptf.SbMax_),
    ScMin_(ptf.ScMin_),
    ScMax_(ptf.ScMax_),
    limitSaturation_(ptf.limitSaturation_),
    capillaryTolerance_(ptf.capillaryTolerance_),
    singularTolerance_(ptf.singularTolerance_)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyNoFluxSaturation::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Saturation fields
    const fvPatchField<scalar>& SaPatch =
        patch().lookupPatchField<volScalarField, scalar>(SaName_);

    const fvPatchField<scalar>& SbPatch =
        patch().lookupPatchField<volScalarField, scalar>(SbName_);

    const scalarField SaP(SaPatch.patchInternalField());
    const scalarField SbP(SbPatch.patchInternalField());


    // Face fields used by the solver
    const fvsPatchField<scalar>& Faf =
        patch().lookupPatchField<surfaceScalarField, scalar>(FafName_);

    const fvsPatchField<scalar>& Fbf =
        patch().lookupPatchField<surfaceScalarField, scalar>(FbfName_);

    const fvsPatchField<scalar>& Laf =
        patch().lookupPatchField<surfaceScalarField, scalar>(LafName_);

    const fvsPatchField<scalar>& Lbf =
        patch().lookupPatchField<surfaceScalarField, scalar>(LbfName_);

    const fvsPatchField<scalar>& Lf =
        patch().lookupPatchField<surfaceScalarField, scalar>(LfName_);

    const fvsPatchField<scalar>& Maf =
        patch().lookupPatchField<surfaceScalarField, scalar>(MafName_);

    const fvsPatchField<scalar>& Mbf =
        patch().lookupPatchField<surfaceScalarField, scalar>(MbfName_);

    const fvsPatchField<scalar>& phiG =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiGName_);


    // Capillary derivatives
    const fvsPatchField<scalar>& dpcgoDsw =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            dpcgoDswName_
        );

    const fvsPatchField<scalar>& dpcgoDsg =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            dpcgoDsgName_
        );

    const fvsPatchField<scalar>& dpcowDsw =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            dpcowDswName_
        );

    const fvsPatchField<scalar>& dpcowDsg =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            dpcowDsgName_
        );


    scalarField& grad = gradient();

    const scalarField& magSf = patch().magSf();
    const scalarField& delta = patch().deltaCoeffs();


    // --------------------------------------------------------------- //
    // Check whether capillary effects are active on this patch.
    //
    // Saturation gradients affect the phase fluxes only through
    // the mobility-weighted capillary derivatives:
    //
    //   phiPc_a = Maf*(dpcgo_dsg*gradSa + dpcgo_dsw*gradSb)*A
    //   phiPc_b = Mbf*(dpcow_dsg*gradSa + dpcow_dsw*gradSb)*A
    //
    // If all these coefficients are negligible on every face,
    // capillarity cannot be used to determine gradSa/gradSb.
    // --------------------------------------------------------------- //

    bool capillaryActive = false;

    forAll(grad, facei)
    {
        if
        (
            mag(Maf[facei]*dpcgoDsg[facei]) > capillaryTolerance_
         || mag(Maf[facei]*dpcgoDsw[facei]) > capillaryTolerance_
         || mag(Mbf[facei]*dpcowDsg[facei]) > capillaryTolerance_
         || mag(Mbf[facei]*dpcowDsw[facei]) > capillaryTolerance_
        )
        {
            capillaryActive = true;
            break;
        }
    }


    // No capillary effects:
    // use a neutral numerical gradient. The physical impermeability
    // condition must then be imposed directly on phia/phib/phic
    // by the solver.
    if (!capillaryActive)
    {
        grad = scalar(0.0);

        fixedGradientFvPatchScalarField::updateCoeffs();

        return;
    }


    forAll(grad, facei)
    {
        const scalar Fa = Faf[facei];
        const scalar Fb = Fbf[facei];

        const scalar Lt = Lf[facei];
        const scalar Ma = Maf[facei];
        const scalar Mb = Mbf[facei];

        if (mag(Lt) <= VSMALL || magSf[facei] <= VSMALL)
        {
            FatalErrorInFunction
                << "Cannot evaluate darcyNoFluxSaturation on patch "
                << patch().name() << ", face " << facei
                << ": Lf or face area is zero."
                << exit(FatalError);
        }

        const scalar alphaA = Laf[facei]/Lt;
        const scalar alphaB = Lbf[facei]/Lt;


        /*
            Solver definitions:

                phi = phiP + phiG + phiPc_a - phiPc_b

                phia =
                    Fa*phiP + alphaA*phiG + phiPc_a

                phib =
                    Fb*phiP + alphaB*phiG - phiPc_b

            The pressure BC imposes phi = 0:

                phiP = -phiG - phiPc_a + phiPc_b

            Substituting into phia = 0:

                (1-Fa)*phiPc_a + Fa*phiPc_b
                  = (Fa-alphaA)*phiG

            Substituting into phib = 0:

                Fb*phiPc_a + (1-Fb)*phiPc_b
                  = (alphaB-Fb)*phiG

            with

                phiPc_a =
                    Ma*(dpcgo_dsg*gradSa
                       +dpcgo_dsw*gradSb)*A

                phiPc_b =
                    Mb*(dpcow_dsg*gradSa
                       +dpcow_dsw*gradSb)*A.
        */

        const scalar C11 =
            (1.0 - Fa)*Ma*dpcgoDsg[facei]
          + Fa*Mb*dpcowDsg[facei];

        const scalar C12 =
            (1.0 - Fa)*Ma*dpcgoDsw[facei]
          + Fa*Mb*dpcowDsw[facei];

        const scalar C21 =
            Fb*Ma*dpcgoDsg[facei]
          + (1.0 - Fb)*Mb*dpcowDsg[facei];

        const scalar C22 =
            Fb*Ma*dpcgoDsw[facei]
          + (1.0 - Fb)*Mb*dpcowDsw[facei];

        const scalar R1 =
            (Fa - alphaA)*phiG[facei]/magSf[facei];

        const scalar R2 =
            (alphaB - Fb)*phiG[facei]/magSf[facei];

        const scalar det =
            C11*C22 - C12*C21;

        const scalar detScale =
            max
            (
                max(mag(C11*C22), mag(C12*C21)),
                VSMALL
            );

        if (mag(det) <= singularTolerance_*detScale)
        {
            // Capillarity may be active on the patch but locally
            // unavailable/degenerate on this face. In that case the
            // saturation gradient cannot uniquely enforce phase no-flow.
            // Use a neutral numerical closure here; the solver must impose
            // the phase fluxes directly on impermeable walls.
            grad[facei] = 0.0;
            continue;
        }

        const scalar gradSaRaw =
            (R1*C22 - C12*R2)/det;

        const scalar gradSbRaw =
            (C11*R2 - R1*C21)/det;


        // Coupled limiter for Sa, Sb and Sc = 1 - Sa - Sb
        scalar theta = 1.0;

        if (limitSaturation_)
        {
            const scalar invDelta = 1.0/delta[facei];

            const scalar dSa = gradSaRaw*invDelta;
            const scalar dSb = gradSbRaw*invDelta;
            const scalar dSc = -(dSa + dSb);

            const scalar Sa0 = SaP[facei];
            const scalar Sb0 = SbP[facei];
            const scalar Sc0 = 1.0 - Sa0 - Sb0;

            if (dSa > SMALL)
            {
                theta = min(theta, (SaMax_ - Sa0)/dSa);
            }
            else if (dSa < -SMALL)
            {
                theta = min(theta, (SaMin_ - Sa0)/dSa);
            }

            if (dSb > SMALL)
            {
                theta = min(theta, (SbMax_ - Sb0)/dSb);
            }
            else if (dSb < -SMALL)
            {
                theta = min(theta, (SbMin_ - Sb0)/dSb);
            }

            if (dSc > SMALL)
            {
                theta = min(theta, (ScMax_ - Sc0)/dSc);
            }
            else if (dSc < -SMALL)
            {
                theta = min(theta, (ScMin_ - Sc0)/dSc);
            }

            theta =
                min
                (
                    max(theta, scalar(0.0)),
                    scalar(1.0)
                );
        }


        const scalar gradSa =
            theta*gradSaRaw;

        const scalar gradSb =
            theta*gradSbRaw;


        if (component_ == "Sa")
        {
            grad[facei] = gradSa;
        }
        else
        {
            grad[facei] = gradSb;
        }
    }


    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::darcyNoFluxSaturation::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);

    os.writeKeyword("component")
        << component_ << token::END_STATEMENT << nl;

    writeEntryIfDifferent<word>(os, "Sa", "Sa", SaName_);
    writeEntryIfDifferent<word>(os, "Sb", "Sb", SbName_);

    writeEntryIfDifferent<word>(os, "Faf", "Faf", FafName_);
    writeEntryIfDifferent<word>(os, "Fbf", "Fbf", FbfName_);

    writeEntryIfDifferent<word>(os, "Laf", "Laf", LafName_);
    writeEntryIfDifferent<word>(os, "Lbf", "Lbf", LbfName_);
    writeEntryIfDifferent<word>(os, "Lf", "Lf", LfName_);

    writeEntryIfDifferent<word>(os, "Maf", "Maf", MafName_);
    writeEntryIfDifferent<word>(os, "Mbf", "Mbf", MbfName_);

    writeEntryIfDifferent<word>(os, "phiG", "phiG", phiGName_);

    writeEntryIfDifferent<word>
    (
        os,
        "dpcgo_dsw",
        "dpcgo_dsw_f",
        dpcgoDswName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "dpcgo_dsg",
        "dpcgo_dsg_f",
        dpcgoDsgName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "dpcow_dsw",
        "dpcow_dsw_f",
        dpcowDswName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "dpcow_dsg",
        "dpcow_dsg_f",
        dpcowDsgName_
    );

    os.writeKeyword("SaMin")
        << SaMin_ << token::END_STATEMENT << nl;

    os.writeKeyword("SaMax")
        << SaMax_ << token::END_STATEMENT << nl;

    os.writeKeyword("SbMin")
        << SbMin_ << token::END_STATEMENT << nl;

    os.writeKeyword("SbMax")
        << SbMax_ << token::END_STATEMENT << nl;

    os.writeKeyword("ScMin")
        << ScMin_ << token::END_STATEMENT << nl;

    os.writeKeyword("ScMax")
        << ScMax_ << token::END_STATEMENT << nl;

    os.writeKeyword("limitSaturation")
        << limitSaturation_ << token::END_STATEMENT << nl;

    os.writeKeyword("capillaryTolerance")
        << capillaryTolerance_ << token::END_STATEMENT << nl;

    os.writeKeyword("singularTolerance")
        << singularTolerance_ << token::END_STATEMENT << nl;

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyNoFluxSaturation
);
}

// ************************************************************************* //
