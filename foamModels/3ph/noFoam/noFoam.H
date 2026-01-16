#ifndef noFoam_H
#define noFoam_H

#include "foamModel.H"

namespace Foam
{

class noFoam
:
    public foamModel
{


public:
    TypeName("noFoam");

    noFoam(const dictionary& dict, foamAuxFields* aux);

    // noFoam functions

    // Correct
    virtual void correct
    (
        volScalarField& kra,
        const volVectorField U,
        const volScalarField& Sa,
        const volScalarField& Sb,
        volScalarField& Sc,
        const surfaceScalarField& phia,
        const volScalarField& eps
    ) const;
};

} // End namespace Foam

#endif
