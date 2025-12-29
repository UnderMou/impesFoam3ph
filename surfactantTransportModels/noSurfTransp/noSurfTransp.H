#ifndef noSurfTransp_H
#define noSurfTransp_H

#include "fvCFD.H"
#include "surfactantTransportModel.H"

namespace Foam
{

class noSurfTransp
:
    public surfactantTransportModel
{

public:
    TypeName("noSurfTransp");

    noSurfTransp(const dictionary& dict, foamAuxFields* aux);

    // noSurfTransp functions

    virtual void correct
    (
        const volScalarField& Sb,
        const surfaceScalarField& phib,
        const volScalarField& eps
    ) const;
};

} // End namespace Foam

#endif
