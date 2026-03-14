#ifndef Mehrabi_H
#define Mehrabi_H

#include "relativePermeabilityModel.H"

namespace Foam
{

class Mehrabi
:
    public relativePermeabilityModel
{
    scalar Sa_min_;
    scalar Sb_min_;
    scalar Sc_min_1_;
    scalar Sc_min_2_;
    scalar a_exp_;
    scalar b_exp_;
    scalar c_exp_;
    scalar kra_max_;
    scalar krb_max_1_;
    scalar krb_max_2_;
    scalar krc_max_;
    scalar Cr_;

public:
    TypeName("Mehrabi");

    Mehrabi(const dictionary& dict);

    scalar correctSor
    (
        scalar Csw
    ) const;

    scalar correctKrw0
    (
        scalar Csw
    ) const;

    virtual void correct
    (
        volScalarField& kra,
        volScalarField& krb,
        volScalarField& krc,
        volScalarField& Sa,
        volScalarField& Sb,
        volScalarField& Csw
    ) const;
};

} // End namespace Foam

#endif
