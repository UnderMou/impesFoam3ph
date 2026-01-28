#ifndef STARS_H
#define STARS_H

#include "foamModel.H"

namespace Foam
{

class STARS
:
    public foamModel
{
    scalar fmmob_;
    scalar SF_;
    scalar sfbet_;
    scalar fmcap_;
    scalar epcap_;
    scalar sigma_ba_;
    scalar fmsurf_;
    scalar epsurf_;
    scalar fmoil_;
    scalar floil_;
    scalar epoil_;

    Switch FshearActive_;
    Switch FsurfActive_;
    Switch FoilActive_;

    word   ncaModel_;

public:
    TypeName("STARS");

    STARS(const dictionary& dict, foamAuxFields* aux);

    // STARS Fi functions
    void correct_Fdry
    (
        volScalarField& Fdry,
        const volScalarField& Sb
    ) const;

    void correct_Nca
    (
        volScalarField& Nca,
        const volVectorField& U,
        const volScalarField& K,
        const volVectorField& gradP
    ) const;

    void correct_Fshear
    (
        volScalarField& Fshear,
        volScalarField& Nca
    ) const;

    void correct_Fsurf
    (
        volScalarField& Fsurf,
        volScalarField& Cs
    ) const;

    void correct_Foil
    (
        volScalarField& Foil,
        volScalarField& Sc
    ) const;

    void correct_MRF
    (
        volScalarField& MRF,
        const volScalarField& Fdry,
        const volScalarField& Fshear,
        const volScalarField& Fsurf,
        const volScalarField& Foil
    ) const;

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
