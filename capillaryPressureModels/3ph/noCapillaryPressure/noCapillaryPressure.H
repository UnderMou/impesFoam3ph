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
        volScalarField& pcwo,
        volScalarField& pcgo,
        volScalarField& dpcwo_dsw,
        volScalarField& dpcwo_dsg,
        volScalarField& dpcgo_dsw,
        volScalarField& dpcgo_dsg,
        const volScalarField& Sa,
        const volScalarField& Sb,
        const volScalarField& zeroField
    ) const;
};

} // End namespace Foam

#endif
