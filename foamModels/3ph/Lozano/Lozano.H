#ifndef Lozano_H
#define Lozano_H

#include "foamModel.H"

namespace Foam
{

class Lozano
:
    public foamModel
{
    scalar MRF_;

public:
    TypeName("Lozano");

    Lozano(const dictionary& dict, foamAuxFields* aux);

    // Correct
    virtual void correct
    (
        volScalarField& kra,
        const volVectorField U,
        const volScalarField& Sa,
        const volScalarField& Sb,
        volScalarField& Sc,
        const surfaceScalarField& phia,
        const volScalarField& eps,
        const volScalarField& K,
        const volVectorField& gradP
    ) const;
};

} // End namespace Foam

#endif
