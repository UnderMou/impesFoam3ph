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
    volScalarField* Foil;
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
    volScalarField* CsEq;
    volScalarField* dCsEqdCs;
    // isotherm
    volScalarField* oneField;

    foamAuxFields()
    :
        MRF(nullptr),
        Nca(nullptr),       // STARS
        Fdry(nullptr),
        Fshear(nullptr),
        Fsurf(nullptr),
        Foil(nullptr),
        nD(nullptr),        // Ashoori
        rgrc(nullptr),
        F(nullptr),
        AcumCoeff_nD(nullptr),
        nDLE(nullptr),
        Cs(nullptr),        // SurfTransp
        Fcsw(nullptr),
        AcumCoeff_csw(nullptr),
        CsEq(nullptr),
        dCsEqdCs(nullptr),
        oneField(nullptr)   // isotherm
    {}

    void printStatus() const
    {
        Info<< "foamAuxFields status:" << nl;
        Info<< "  MRF: " << (MRF ? "initialized" : "null") << nl;
        Info<< "  Nca: " << (Nca ? "initialized" : "null") << nl;
        Info<< "  Fdry: " << (Fdry ? "initialized" : "null") << nl;
        Info<< "  Fshear: " << (Fshear ? "initialized" : "null") << nl;
        Info<< "  Fsurf: " << (Fsurf ? "initialized" : "null") << nl;
        Info<< "  Foil: " << (Foil ? "initialized" : "null") << nl;
        Info<< "  nD: " << (nD ? "initialized" : "null") << nl;
        Info<< "  rgrc: " << (rgrc ? "initialized" : "null") << nl;
        Info<< "  F: " << (F ? "initialized" : "null") << nl;
        Info<< "  AcumCoeff_nD: " << (AcumCoeff_nD ? "initialized" : "null") << nl;
        Info<< "  nDLE: " << (nDLE ? "initialized" : "null") << nl;
        Info<< "  Cs: " << (Cs ? "initialized" : "null") << nl;
        Info<< "  Fcsw: " << (Fcsw ? "initialized" : "null") << nl;
        Info<< "  AcumCoeff_csw: " << (AcumCoeff_csw ? "initialized" : "null") << nl;
        Info<< "  CsEq: " << (CsEq ? "initialized" : "null") << nl;
        Info<< "  dCsEqdCs: " << (dCsEqdCs ? "initialized" : "null") << nl;
        Info<< "  oneField: " << (oneField ? "initialized" : "null") << nl;
    }
};

} // End namespace Foam

#endif
