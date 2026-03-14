#ifndef noAdsorption_H
#define noAdsorption_H

#include "isothermModel.H"

namespace Foam
{

class noAdsorption
:
    public isothermModel
{

public:
    TypeName("noAdsorption");

    noAdsorption(const dictionary& dict, foamAuxFields* aux);

    virtual void correct
    (
       
    ) const;
};

} // End namespace Foam

#endif
