#ifndef Langmuir_H
#define Langmuir_H

#include "isothermModel.H"

namespace Foam
{

class Langmuir
:
    public isothermModel
{
    scalar qo_;
    scalar KL_;

public:
    TypeName("Langmuir");

    Langmuir(const dictionary& dict, foamAuxFields* aux);

    virtual void correct
    (
       
    ) const;
};

} // End namespace Foam

#endif
