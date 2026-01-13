#ifndef ConstantInjection_H
#define ConstantInjection_H

#include "wellModel.H"

namespace Foam
{

class ConstantInjection
:
    public wellModel
{
    scalar Fb_inj_;

public:
    TypeName("ConstantInjection");

    ConstantInjection(const dictionary& dict);

    virtual void correct
    (
        volScalarField& qb,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod
    ) const;
};

} // End namespace Foam

#endif
