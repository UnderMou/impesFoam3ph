#ifndef relativePermeabilityModel_H
#define relativePermeabilityModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class relativePermeabilityModel
{
protected:
    const dictionary& dict_;

public:
    TypeName("relativePermeabilityModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        relativePermeabilityModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    relativePermeabilityModel(const dictionary& dict)
    :
        dict_(dict)
    {}

    virtual ~relativePermeabilityModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        volScalarField& kra,
        volScalarField& krb,
        volScalarField& krc,
        volScalarField& Sa,
        volScalarField& Sb
    ) const = 0;

    // Factory
    static autoPtr<relativePermeabilityModel> New
    (
        const dictionary& transportProperties
    );
};

} // End namespace Foam

#endif
