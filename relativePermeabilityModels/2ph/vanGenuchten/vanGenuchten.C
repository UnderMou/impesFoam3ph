#include "vanGenuchten.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    vanGenuchten,
    0
);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    vanGenuchten,
    dictionary
);

vanGenuchten::vanGenuchten(const dictionary& dict)
:
    relativePermeabilityModel(dict),
    Sa_min_(readScalar(dict.lookup("Sa_min"))),
    Sb_min_(readScalar(dict.lookup("Sb_min"))),
    n_(readScalar(dict.lookup("n"))),
    m_(1.0 - 1.0/n_),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_(readScalar(dict.lookup("krb_max")))
{}

void vanGenuchten::correct
(
    volScalarField& kra,
    volScalarField& krb,
    volScalarField& Sb
) const
{

    volScalarField Se = (Sb - Sb_min_) / (1.0 - Sa_min_ - Sb_min_);

    kra = kra_max_ * Foam::pow(scalar(1.0) - Se, 0.5) * Foam::pow(scalar(1.0) - Foam::pow(Se, 1.0/m_), 2*m_);
    krb = krb_max_ * Foam::pow(Se, 0.5) * Foam::pow(scalar(1.0) - Foam::pow(scalar(1.0) - Foam::pow(Se, 1.0/m_), m_), 2); 

}

} // End namespace Foam
