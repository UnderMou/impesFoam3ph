#ifndef Henry_H
#define Henry_H

#include "isothermModel.H"

namespace Foam
{

class Henry
:
    public isothermModel
{
    scalar Kd_;

public:
    TypeName("Henry");

    Henry(const dictionary& dict, foamAuxFields* aux);

    virtual void correct
    (
       
    ) const;
};

} // End namespace Foam

#endif
