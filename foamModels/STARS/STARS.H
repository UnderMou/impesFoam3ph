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

public:
    TypeName("STARS");

    STARS(const dictionary& dict);

    // STARS Fi functions
    void correct_Fdry
    (
        volScalarField& Fdry,
        const volScalarField Sb
    ) const;

    void correct_Nca
    (
        volScalarField& Nca,
        const volVectorField U,
        const dimensionedScalar mu_b
    ) const;

    void correct_Fshear
    (
        volScalarField& Fshear,
        volScalarField& Nca
    ) const;

    void correct_MRF
    (
        volScalarField& MRF,
        const volScalarField Fdry,
        const volScalarField Fshear
    ) const;

    // Correct
    virtual void correct
    (
        volScalarField& kra,
        volScalarField& MRF,
        volScalarField& Fdry,
        volScalarField& Fshear,
        volScalarField& Nca,
        const volVectorField U,
        const volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
