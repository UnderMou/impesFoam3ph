#ifndef foamAuxFields_H
#define foamAuxFields_H

#include "volFields.H"

namespace Foam
{

struct foamAuxFields
{
    volScalarField* MRF;
    // STARS
    volScalarField* Nca;
    volScalarField* Fdry;
    volScalarField* Fshear;
    // Ashoori
    volScalarField* nD;
    volScalarField* rgrc;
    volScalarField* F;
    volScalarField* AcumCoeff;
    volScalarField* nDLE;

    foamAuxFields()
    :
        MRF(nullptr),
        Nca(nullptr),   // STARS
        Fdry(nullptr),
        Fshear(nullptr),
        nD(nullptr),    // Ashoori
        rgrc(nullptr),
        F(nullptr),
        AcumCoeff(nullptr),
        nDLE(nullptr)
    {}
};

} // End namespace Foam

#endif
