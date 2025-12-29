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
    volScalarField* Fsurf;
    // Ashoori
    volScalarField* nD;
    volScalarField* rgrc;
    volScalarField* F;
    volScalarField* AcumCoeff_nD;
    volScalarField* nDLE;
    // SurfTransp
    volScalarField* Cs;
    volScalarField* Fcsw;
    volScalarField* AcumCoeff_csw;

    foamAuxFields()
    :
        MRF(nullptr),
        Nca(nullptr),   // STARS
        Fdry(nullptr),
        Fshear(nullptr),
        Fsurf(nullptr),
        nD(nullptr),    // Ashoori
        rgrc(nullptr),
        F(nullptr),
        AcumCoeff_nD(nullptr),
        nDLE(nullptr),
        Cs(nullptr),    // SurfTransp
        Fcsw(nullptr),
        AcumCoeff_csw(nullptr)
    {}
};

} // End namespace Foam

#endif
