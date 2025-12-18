#ifndef foamModel_H
#define foamModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class foamModel
{
protected:
    const dictionary& dict_;
    dimensionedScalar mu_b_;

public:
    TypeName("foamModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        foamModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    foamModel(const dictionary& dict)
    :
        dict_(dict),
        mu_b_(dict.parent().lookup("mu_b"))
    {}

    virtual ~foamModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        volScalarField& kra,
        volScalarField& MRF,
        volScalarField& Fdry,
        volScalarField& Fshear,
        volScalarField& Nca,
        const volVectorField U,
        const volScalarField& Sb
    ) const = 0;

    // Factory
    static autoPtr<foamModel> New
    (
        const dictionary& transportProperties
    );
};

} // End namespace Foam

#endif
