#include "ashoori.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    ashoori,
    0
);

addToRunTimeSelectionTable
(
    foamModel,
    ashoori,
    dictionary
);

ashoori::ashoori(const dictionary& dict, foamAuxFields* aux)
:
    foamModel(dict, aux),
    Cmrf_(readScalar(dict.lookup("Cmrf"))),
    Kc_(readScalar(dict.lookup("Kc"))),
    nMAX_(readScalar(dict.lookup("nMAX"))),
    Sw_ast_(readScalar(dict.lookup("Sw_ast"))),
    A_(readScalar(dict.lookup("A")))
{}

void ashoori::correct_nD_LE
(
    volScalarField& nDLE,
    volScalarField Sb
) const
{
    // Correct at cell center
    forAll(nDLE,i){
        if (Sb[i] > Sw_ast_){
            nDLE[i] = Foam::tanh(A_ * (Sb[i] - Sw_ast_));
        }
    }

    // Correct at the boundaries
    forAll(Sb.boundaryField(), patchi)
    {
        fvPatchScalarField& nDLEPatch = nDLE.boundaryFieldRef()[patchi];
        fvPatchScalarField& SbPatch = Sb.boundaryFieldRef()[patchi];

        forAll(SbPatch, facei)
        {
            if (SbPatch[facei] > Sw_ast_){
                nDLEPatch[facei] = Foam::tanh(A_ * (SbPatch[facei] - Sw_ast_));
            }
        }
    }
}

void ashoori::correct_rgrc
(
    volScalarField& rgrc,
    volScalarField& nDLE,
    const volScalarField Sb,
    const volScalarField nD
) const
{
    correct_nD_LE(nDLE, Sb);
    rgrc = Kc_ * nMAX_ * (nDLE - nD);
}

void ashoori::correct_MRF
(
    volScalarField& MRF,
    const volScalarField nD
) const
{
    MRF = Cmrf_ * nD + scalar(1.0);
}

void ashoori::correct
(
    volScalarField& kra,
    const volVectorField U,
    const volScalarField& Sa,
    const volScalarField& Sb,
    const surfaceScalarField& phia,
    const volScalarField& eps
) const
{
    // Info << "fmmob = " << fmmob_ << endl;
    // Info << "mu_w = " << mu_b_.value() << endl;

    if (!aux_ || !aux_->rgrc)
    {
        FatalErrorInFunction
            << "ashoori requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& nD        = *aux_->nD;
    volScalarField& rgrc      = *aux_->rgrc;
    volScalarField& F         = *aux_->F;
    volScalarField& AcumCoeff = *aux_->AcumCoeff_nD;
    volScalarField& nDLE      = *aux_->nDLE;
    volScalarField& MRF       = *aux_->MRF;

    F = -fvc::div(phia); // + qa
    AcumCoeff = eps*(Sa + VSMALL);
    correct_rgrc(rgrc, nDLE, Sb, nD);

    Info << "rgrc OK" << endl;
    // printField(rgrc);
    
    fvScalarMatrix nDEqn
    (
        AcumCoeff*fvm::ddt(nD) + fvc::div(phia, nD) 
        ==
        (eps*(Sa + VSMALL)/nMAX_)*rgrc + fvc::Sp(-F, nD) // + qf
    );

    nDEqn.solve();

    correct_MRF(MRF, nD);
    
    kra /= MRF;
}

} // End namespace Foam
