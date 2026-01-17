#ifndef wellModel_H
#define wellModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class wellModel
{
protected:
    const dictionary& dict_;

public:
    TypeName("wellModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        wellModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    wellModel(const dictionary& dict)
    :
        dict_(dict)
    {}

    virtual ~wellModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        volScalarField& qb,
        const volScalarField& Fb, 
        const volScalarField& qt_inj,
        const volScalarField& qt_prod
    ) const = 0;

    // Factory
    static autoPtr<wellModel> New
    (
        const dictionary& wellProperties
    );
};

} // End namespace Foam

#endif
