#ifndef Csw_H
#define Csw_H

#include "fvCFD.H"
#include "surfactantTransportModel.H"

namespace Foam
{

class Csw
:
    public surfactantTransportModel
{
    scalar Kd_;
    scalar rho_s_;

    scalar rho_sw_;

public:
    TypeName("Csw");

    Csw(const dictionary& dict, foamAuxFields* aux);

    // Csw functions

    scalar correct_dCsEqdCs
    (

    ) const;

    virtual void correct
    (
        const volScalarField& Sb,
        const surfaceScalarField& phib,
        const volScalarField& eps
    ) const;
};

} // End namespace Foam

#endif
