#ifndef isothermModel_H
#define isothermModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "foamAuxFields.H"

namespace Foam
{

class isothermModel
{
protected:
    const dictionary& dict_;
    foamAuxFields* aux_;

public:
    TypeName("isothermModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        isothermModel,
        dictionary,
        (
            const dictionary& dict,
            foamAuxFields* aux
        ),
        (dict, aux)
    );

    isothermModel(const dictionary& dict, foamAuxFields* aux)
    :
        dict_(dict),
        aux_(aux)
    {}

    virtual ~isothermModel() = default;

    // ---- Interface ----

    virtual void correct
    (

    ) const = 0;

    // Factory
    static autoPtr<isothermModel> New
    (
        const dictionary& transportProperties,
        foamAuxFields* aux
    );
};

} // End namespace Foam

#endif
