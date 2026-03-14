#include "Mehrabi.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    Mehrabi,
    0
);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    Mehrabi,
    dictionary
);

Mehrabi::Mehrabi(const dictionary& dict)
:
    relativePermeabilityModel(dict),
    Sa_min_(readScalar(dict.lookup("Sa_min"))),
    Sb_min_(readScalar(dict.lookup("Sb_min"))),
    Sc_min_1_(readScalar(dict.lookup("Sc_min_1"))),
    Sc_min_2_(readScalar(dict.lookup("Sc_min_2"))),
    a_exp_(readScalar(dict.lookup("a_exp"))),
    b_exp_(readScalar(dict.lookup("b_exp"))),
    c_exp_(readScalar(dict.lookup("c_exp"))),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_1_(readScalar(dict.lookup("krb_max_1"))),
    krb_max_2_(readScalar(dict.lookup("krb_max_2"))),
    krc_max_(readScalar(dict.lookup("krc_max"))),
    Cr_(readScalar(dict.lookup("Cr")))
{}

scalar Mehrabi::correctSor
(
    scalar Csw
) const
{
    if (Csw <= Cr_){
        return Sc_min_2_;
    }
    else{
        return Sc_min_1_;
    }
}

scalar Mehrabi::correctKrw0
(
    scalar Csw
) const
{
    if (Csw <= Cr_){
        return krb_max_2_;
    }
    else{
        return krb_max_1_;
    }
}

void Mehrabi::correct
(
    volScalarField& kra,
    volScalarField& krb,
    volScalarField& krc,
    volScalarField& Sa,
    volScalarField& Sb,
    volScalarField& Csw
) const
{

    volScalarField Sc = scalar(1.0) - Sa - Sb;

    scalar Sc_min_corr = Sc_min_2_;
    scalar krb_max_corr = krb_max_2_;

    scalar Se_a = (Sa[0] - Sa_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
    scalar Se_b = (Sb[0] - Sb_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
    scalar Se_c = (Sc[0] - Sc_min_corr) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);

    // Correct kra
    // cell centered values
    forAll(kra,i){

        Sc_min_corr = correctSor(Csw[i]);

        if (Sa[i]<=Sa_min_ + SMALL){
            kra[i]=0.0;
        }
        else if (Sa[i]>=1.0 - Sb_min_ - Sc_min_corr - SMALL){
            kra[i]=kra_max_;
        }
        else{
            Se_a = (Sa[i] - Sa_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
            kra[i] = kra_max_ * Foam::pow(Se_a, a_exp_);
        }
    }
    // face centered values
    forAll(kra.boundaryField(), patchi)
    {
        fvPatchScalarField& kraPatch = kra.boundaryFieldRef()[patchi];
        fvPatchScalarField& SaPatch = Sa.boundaryFieldRef()[patchi];
        fvPatchScalarField& CswPatch = Csw.boundaryFieldRef()[patchi];

        forAll(kraPatch, facei)
        {
            Sc_min_corr = correctSor(CswPatch[facei]);

            if (SaPatch[facei]<=Sa_min_ + SMALL){
                kraPatch[facei]=0.0;
            }
            else if (SaPatch[facei]>=1.0 - Sb_min_ - Sc_min_corr - SMALL){
                kraPatch[facei]=kra_max_;
            }
            else{
                Se_a = (SaPatch[facei] - Sa_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
                kraPatch[facei] = kra_max_ * Foam::pow(Se_a, a_exp_);
            }
        }
    }

    // Correct krb
    // cell centered values
    forAll(krb,i){

        Sc_min_corr = correctSor(Csw[i]);
        krb_max_corr = correctKrw0(Csw[i]);

        if (Sb[i]<=Sb_min_ + SMALL){
            krb[i]=0.0;
        }
        else if (Sb[i]>=1.0 - Sa_min_ - Sc_min_corr - SMALL){
            krb[i]=krb_max_corr;
        }
        else{
            Se_b = (Sb[i] - Sb_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
            krb[i] = krb_max_corr * Foam::pow(Se_b, b_exp_);
        }
    }
    // face centered values
    forAll(krb.boundaryField(), patchi)
    {
        fvPatchScalarField& krbPatch = krb.boundaryFieldRef()[patchi];
        fvPatchScalarField& SbPatch = Sb.boundaryFieldRef()[patchi];
        fvPatchScalarField& CswPatch = Csw.boundaryFieldRef()[patchi];

        forAll(krbPatch, facei)
        {
            Sc_min_corr = correctSor(CswPatch[facei]);
            krb_max_corr = correctKrw0(CswPatch[facei]);

            if (SbPatch[facei]<=Sb_min_ + SMALL){
                krbPatch[facei]=0.0;
            }
            else if (SbPatch[facei]>=1.0 - Sa_min_ - Sc_min_corr - SMALL){
                krbPatch[facei]=krb_max_corr;
            }
            else{
                Se_b = (SbPatch[facei] - Sb_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
                krbPatch[facei] = krb_max_corr * Foam::pow(Se_b, b_exp_);
            }
        }
    }

    // Correct krc
    // cell centered values
    forAll(krc,i){

        Sc_min_corr = correctSor(Csw[i]);

        if (Sc[i]<=Sc_min_corr + SMALL){
            krc[i]=0.0;
        }
        else if (Sc[i]>=1.0 - Sb_min_ - Sa_min_ - SMALL){
            krc[i]=krc_max_;
        }
        else{
            Se_c = (Sc[i] - Sc_min_corr) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
            krc[i] = krc_max_ * Foam::pow(Se_c, c_exp_);
        }
    }
    // face centered values
    forAll(krc.boundaryField(), patchi)
    {
        fvPatchScalarField& krcPatch = krc.boundaryFieldRef()[patchi];
        fvPatchScalarField& ScPatch = Sc.boundaryFieldRef()[patchi];
        fvPatchScalarField& CswPatch = Csw.boundaryFieldRef()[patchi];

        forAll(krcPatch, facei)
        {
            Sc_min_corr = correctSor(CswPatch[facei]);

            if (ScPatch[facei]<=Sc_min_corr + SMALL){
                krcPatch[facei]=0.0;
            }
            else if (ScPatch[facei]>=1.0 - Sb_min_ - Sa_min_ - SMALL){
                krcPatch[facei]=krc_max_;
            }
            else{
                Se_c = (ScPatch[facei] - Sc_min_corr) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_corr);
                krcPatch[facei] = krc_max_ * Foam::pow(Se_c, c_exp_);
            }
        }
    }
}

} // End namespace Foam
