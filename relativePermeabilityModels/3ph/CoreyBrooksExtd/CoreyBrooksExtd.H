#ifndef CoreyBrooksExtd_H
#define CoreyBrooksExtd_H

#include "relativePermeabilityModel.H"

namespace Foam
{

class CoreyBrooksExtd
:
    public relativePermeabilityModel
{
    scalar Sa_min_;
    scalar Sb_min_;
    scalar Sc_min_;
    scalar a_exp_;
    scalar b_exp_;
    scalar c_exp_;
    scalar kra_max_;
    scalar krb_max_;
    scalar krc_max_;

public:
    TypeName("CoreyBrooksExtd");

    CoreyBrooksExtd(const dictionary& dict);

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
