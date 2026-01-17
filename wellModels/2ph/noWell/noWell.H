#ifndef noWell_H
#define noWell_H

#include "wellModel.H"

namespace Foam
{

class noWell
:
    public wellModel
{

public:
    TypeName("noWell");

    noWell(const dictionary& dict);

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
