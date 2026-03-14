#include "CoreyBrooks.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    CoreyBrooks,
    0
);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    CoreyBrooks,
    dictionary
);

CoreyBrooks::CoreyBrooks(const dictionary& dict)
:
    relativePermeabilityModel(dict),
    Sa_min_(readScalar(dict.lookup("Sa_min"))),
    Sb_min_(readScalar(dict.lookup("Sb_min"))),
    Sc_min_(readScalar(dict.lookup("Sc_min"))),
    a_exp_(readScalar(dict.lookup("a_exp"))),
    b_exp_(readScalar(dict.lookup("b_exp"))),
    c_exp_(readScalar(dict.lookup("c_exp"))),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_(readScalar(dict.lookup("krb_max"))),
    krc_max_(readScalar(dict.lookup("krc_max")))
{}

void CoreyBrooks::correct
(
    volScalarField& kra,
    volScalarField& krb,
    volScalarField& krc,
    volScalarField& Sa,
    volScalarField& Sb,
    volScalarField& Csw
) const
{

    volScalarField Se_a = (Sa - Sa_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_);
    volScalarField Se_b = (Sb - Sb_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_);
    volScalarField Se_c = (scalar(1.0) - Sa - Sb - Sc_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_);
    
    kra = kra_max_ * Foam::pow(Se_a, a_exp_); 
    krb = krb_max_ * Foam::pow(Se_b, b_exp_); 
    krc = krc_max_ * Foam::pow(Se_c, c_exp_); 

}

} // End namespace Foam
