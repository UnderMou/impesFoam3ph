#ifndef CoreyBrooksPC_H
#define CoreyBrooksPC_H

#include "capillaryPressureModel.H"

namespace Foam
{

class CoreyBrooksPC
:
    public capillaryPressureModel
{
    scalar epsilon_;

public:
    TypeName("CoreyBrooksPC");

    CoreyBrooksPC(const dictionary& dict);

    virtual void correct
    (
        volScalarField& pcow,
        volScalarField& pcgo,
        volScalarField& dpcow_dsw,
        volScalarField& dpcow_dsg,
        volScalarField& dpcgo_dsw,
        volScalarField& dpcgo_dsg,
        const volScalarField& Sa,
        const volScalarField& Sb,
        const volScalarField& zeroField
    ) const;
};

} // End namespace Foam

#endif
