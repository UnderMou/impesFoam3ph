#ifndef MehrabiFoam_H
#define MehrabiFoam_H

#include "foamModel.H"

namespace Foam
{

class MehrabiFoam
:
    public foamModel
{
    scalar MRF_;
    scalar Cr_;

public:
    TypeName("MehrabiFoam");

    MehrabiFoam(const dictionary& dict, foamAuxFields* aux);

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
