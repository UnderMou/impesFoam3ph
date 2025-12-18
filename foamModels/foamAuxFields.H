#ifndef foamAuxFields_H
#define foamAuxFields_H

#include "volFields.H"

namespace Foam
{

struct foamAuxFields
{
    volScalarField* Nca;
    volScalarField* Fdry;
    volScalarField* Fshear;
    volScalarField* MRF;

    foamAuxFields()
    :
        Nca(nullptr),
        Fdry(nullptr),
        Fshear(nullptr),
        MRF(nullptr)
    {}
};

} // End namespace Foam

#endif
