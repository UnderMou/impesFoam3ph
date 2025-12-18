#ifndef CoreyBrooks_H
#define CoreyBrooks_H

#include "relativePermeabilityModel.H"

namespace Foam
{

class CoreyBrooks
:
    public relativePermeabilityModel
{
    scalar Sa_min_;
    scalar Sb_min_;
    scalar a_exp_;
    scalar b_exp_;
    scalar kra_max_;
    scalar krb_max_;

public:
    TypeName("CoreyBrooks");

    CoreyBrooks(const dictionary& dict);

    virtual void correct
    (
        volScalarField& kra,
        volScalarField& krb,
        const volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
