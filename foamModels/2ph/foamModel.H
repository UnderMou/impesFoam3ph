#ifndef foamModel_H
#define foamModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "foamAuxFields.H"

namespace Foam
{

class foamModel
{
protected:
    const dictionary& dict_;
    foamAuxFields* aux_;

    dimensionedScalar mu_b_;

public:
    TypeName("foamModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        foamModel,
        dictionary,
        (
            const dictionary& dict,
            foamAuxFields* aux
        ),
        (dict, aux)
    );

    foamModel(const dictionary& dict, foamAuxFields* aux)
    :
        dict_(dict),
        aux_(aux),
        mu_b_(dict.parent().lookup("mu_b"))
    {}

    virtual ~foamModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        volScalarField& kra,
        const volVectorField U,
        const volScalarField& Sa,
        const volScalarField& Sb,
        const surfaceScalarField& phia,
        const volScalarField& eps
    ) const = 0;

    // Factory
    static autoPtr<foamModel> New
    (
        const dictionary& transportProperties,
        foamAuxFields* aux
    );
};

} // End namespace Foam

#endif
