#ifndef vanGenuchten_H
#define vanGenuchten_H

#include "relativePermeabilityModel.H"

namespace Foam
{

class vanGenuchten
:
    public relativePermeabilityModel
{
    scalar Sa_min_;
    scalar Sb_min_;
    scalar n_;
    scalar m_;
    scalar kra_max_;
    scalar krb_max_;

public:
    TypeName("vanGenuchten");

    vanGenuchten(const dictionary& dict);

    virtual void correct
    (
        volScalarField& kra,
        volScalarField& krb,
        volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
