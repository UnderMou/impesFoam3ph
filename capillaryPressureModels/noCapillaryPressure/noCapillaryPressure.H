#ifndef noCapillaryPressure_H
#define noCapillaryPressure_H

#include "capillaryPressureModel.H"

namespace Foam
{

class noCapillaryPressure
:
    public capillaryPressureModel
{

public:
    TypeName("noCapillaryPressure");

    noCapillaryPressure(const dictionary& dict);

    virtual void correct
    (
        volScalarField& pc,
        volScalarField& dpcds,
        const volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
