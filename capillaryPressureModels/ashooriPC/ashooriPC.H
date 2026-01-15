#ifndef ashooriPC_H
#define ashooriPC_H

#include "capillaryPressureModel.H"

namespace Foam
{

class ashooriPC
:
    public capillaryPressureModel
{
    scalar Swc_;
    scalar Sgr_;
    scalar pc0_;
    scalar c_;

public:
    TypeName("ashooriPC");

    ashooriPC(const dictionary& dict);

    virtual void correct
    (
        volScalarField& pc,
        volScalarField& dpcds,
        const volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
