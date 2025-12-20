#ifndef CoreyBrooksPC_H
#define CoreyBrooksPC_H

#include "capillaryPressureModel.H"

namespace Foam
{

class CoreyBrooksPC
:
    public capillaryPressureModel
{
    scalar Spc_irr_;
    scalar Spc_max_;
    scalar pc0_;
    double alpha_;

public:
    TypeName("CoreyBrooksPC");

    CoreyBrooksPC(const dictionary& dict);

    void correct_dpcds
    (
        volScalarField& dpcds,
        const volScalarField& Sb
    ) const;

    virtual void correct
    (
        volScalarField& pc,
        volScalarField& dpcds,
        const volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
