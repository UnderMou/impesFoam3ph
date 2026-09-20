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
    McfName_("Mcf"),
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
    mobilityTolerance_(1e-12),
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
    McfName_(dict.lookupOrDefault<word>("Mcf", "Mcf")),
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
    mobilityTolerance_
    (
        dict.lookupOrDefault<scalar>("mobilityTolerance", 1e-12)
    ),
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
    McfName_(ptf.McfName_),
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
    mobilityTolerance_(ptf.mobilityTolerance_),
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
    McfName_(ptf.McfName_),
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
    mobilityTolerance_(ptf.mobilityTolerance_),
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

    const fvsPatchField<scalar>& Mcf =
        patch().lookupPatchField<surfaceScalarField, scalar>(McfName_);

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
        const scalar Mc = Mcf[facei];

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

            Hence the two phase no-flow relations are

                (1-Fa)*phiPc_a + Fa*phiPc_b
                    = (Fa-alphaA)*phiG                          (A)

                Fb*phiPc_a + (1-Fb)*phiPc_b
                    = (alphaB-Fb)*phiG                         (B)

            with

                phiPc_a =
                    Ma*(dpcgo_dsg*gradSa
                       +dpcgo_dsw*gradSb)*A

                phiPc_b =
                    Mb*(dpcow_dsg*gradSa
                       +dpcow_dsw*gradSb)*A.

            When all three phases are mobile, (A) and (B) form a 2x2
            system for gradSa and gradSb.

            If one phase becomes immobile, only one independent no-flow
            relation remains. In that case a minimum-norm solution of the
            remaining scalar equation is used instead of setting both
            saturation gradients to zero.
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


        // ----------------------------------------------------------- //
        // Determine which phases are locally mobile.
        //
        // mobilityTolerance_ is RELATIVE to the largest phase mobility
        // on this face, avoiding an absolute dimensional mobility cutoff.
        // ----------------------------------------------------------- //

        const scalar mobilityScale =
            max
            (
                max(mag(Ma), mag(Mb)),
                mag(Mc)
            );

        if (mobilityScale <= VSMALL)
        {
            grad[facei] = 0.0;
            continue;
        }

        const scalar mobilityCutoff =
            mobilityTolerance_*mobilityScale;

        const bool aActive = mag(Ma) > mobilityCutoff;
        const bool bActive = mag(Mb) > mobilityCutoff;
        const bool cActive = mag(Mc) > mobilityCutoff;

        const label nActive =
            label(aActive) + label(bActive) + label(cActive);


        scalar gradSaRaw = 0.0;
        scalar gradSbRaw = 0.0;
        bool gradientDetermined = false;


        // ----------------------------------------------------------- //
        // Three mobile phases: solve the full 2x2 system.
        // ----------------------------------------------------------- //

        if (nActive == 3)
        {
            const scalar det =
                C11*C22 - C12*C21;

            const scalar detScale =
                max
                (
                    max(mag(C11*C22), mag(C12*C21)),
                    VSMALL
                );

            if (mag(det) > singularTolerance_*detScale)
            {
                gradSaRaw =
                    (R1*C22 - C12*R2)/det;

                gradSbRaw =
                    (C11*R2 - R1*C21)/det;

                gradientDetermined = true;
            }
            else
            {
                // Rank-deficient capillary response even though all phases
                // are mobile. Enforce the better-conditioned independent
                // relation using its minimum-norm solution.
                const scalar rowNormA2 = C11*C11 + C12*C12;
                const scalar rowNormB2 = C21*C21 + C22*C22;

                if
                (
                    rowNormA2 >= rowNormB2
                 && rowNormA2 > VSMALL
                )
                {
                    gradSaRaw = R1*C11/rowNormA2;
                    gradSbRaw = R1*C12/rowNormA2;
                    gradientDetermined = true;
                }
                else if (rowNormB2 > VSMALL)
                {
                    gradSaRaw = R2*C21/rowNormB2;
                    gradSbRaw = R2*C22/rowNormB2;
                    gradientDetermined = true;
                }
            }
        }


        // ----------------------------------------------------------- //
        // Phase a immobile: phia is already zero to mobility accuracy.
        // Enforce the remaining b-c no-flow relation (B).
        // ----------------------------------------------------------- //

        else if (!aActive && bActive && cActive)
        {
            const scalar rowNorm2 =
                C21*C21 + C22*C22;

            if (rowNorm2 > VSMALL)
            {
                gradSaRaw = R2*C21/rowNorm2;
                gradSbRaw = R2*C22/rowNorm2;
                gradientDetermined = true;
            }
        }


        // ----------------------------------------------------------- //
        // Phase b immobile: phib is already zero to mobility accuracy.
        // Enforce the remaining a-c no-flow relation (A).
        // ----------------------------------------------------------- //

        else if (aActive && !bActive && cActive)
        {
            const scalar rowNorm2 =
                C11*C11 + C12*C12;

            if (rowNorm2 > VSMALL)
            {
                gradSaRaw = R1*C11/rowNorm2;
                gradSbRaw = R1*C12/rowNorm2;
                gradientDetermined = true;
            }
        }


        // ----------------------------------------------------------- //
        // Phase c immobile: only phases a and b are mobile.
        //
        // Since phi = phia + phib = 0, imposing either one of the two
        // phase no-flow equations automatically imposes the other.
        // Choose the better-conditioned row and use its minimum-norm
        // solution.
        // ----------------------------------------------------------- //

        else if (aActive && bActive && !cActive)
        {
            const scalar rowNormA2 =
                C11*C11 + C12*C12;

            const scalar rowNormB2 =
                C21*C21 + C22*C22;

            if
            (
                rowNormA2 >= rowNormB2
             && rowNormA2 > VSMALL
            )
            {
                gradSaRaw = R1*C11/rowNormA2;
                gradSbRaw = R1*C12/rowNormA2;
                gradientDetermined = true;
            }
            else if (rowNormB2 > VSMALL)
            {
                gradSaRaw = R2*C21/rowNormB2;
                gradSbRaw = R2*C22/rowNormB2;
                gradientDetermined = true;
            }
        }


        // ----------------------------------------------------------- //
        // Zero or one mobile phase:
        //
        // The pressure condition phi = 0 already guarantees no-flow of
        // the only mobile phase. Impermeability does not determine Sa/Sb
        // gradients, so use a neutral numerical closure.
        // ----------------------------------------------------------- //

        if (!gradientDetermined)
        {
            grad[facei] = 0.0;
            continue;
        }


        // ----------------------------------------------------------- //
        // Coupled limiter for Sa, Sb and Sc = 1 - Sa - Sb.
        // Both gradients are multiplied by the SAME theta.
        // ----------------------------------------------------------- //

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
    writeEntryIfDifferent<word>(os, "Mcf", "Mcf", McfName_);

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

    os.writeKeyword("mobilityTolerance")
        << mobilityTolerance_ << token::END_STATEMENT << nl;

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
