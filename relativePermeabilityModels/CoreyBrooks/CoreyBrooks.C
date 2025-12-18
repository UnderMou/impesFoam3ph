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
    a_exp_(readScalar(dict.lookup("a_exp"))),
    b_exp_(readScalar(dict.lookup("b_exp"))),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_(readScalar(dict.lookup("krb_max")))
{}

void CoreyBrooks::correct
(
    volScalarField& kra,
    volScalarField& krb,
    const volScalarField& Sb
) const
{

    volScalarField Se = (Sb - Sb_min_) / (1.0 - Sa_min_ - Sb_min_);

    kra = kra_max_ * Foam::pow(scalar(1.0) - Se, a_exp_);
    krb = krb_max_ * Foam::pow(Se, b_exp_); 

}

} // End namespace Foam
