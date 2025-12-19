#ifndef ashoori_H
#define ashoori_H

#include "fvCFD.H"
#include "foamModel.H"

namespace Foam
{

class ashoori
:
    public foamModel
{
    scalar Cmrf_;
    scalar Kc_;
    scalar nMAX_;
    scalar Sw_ast_;
    scalar A_;

public:
    TypeName("ashoori");

    ashoori(const dictionary& dict, foamAuxFields* aux);

    // ashoori functions
    void correct_nD_LE
    (
        volScalarField& nDLE,
        volScalarField Sb
    ) const;

    void correct_rgrc
    (
        volScalarField& rgrc,
        volScalarField& nDLE,
        const volScalarField Sb,
        const volScalarField nD
    ) const;

    void correct_MRF
    (
        volScalarField& MRF,
        const volScalarField nD
    ) const;

    // Correct
    virtual void correct
    (
        volScalarField& kra,
        const volVectorField U,
        const volScalarField& Sa,
        const volScalarField& Sb,
        const surfaceScalarField& phia,
        const volScalarField& eps
    ) const;
};

} // End namespace Foam

#endif
