#ifndef capillaryPressureModel_H
#define capillaryPressureModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class capillaryPressureModel
{
protected:
    const dictionary& dict_;

public:
    TypeName("capillaryPressureModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        capillaryPressureModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    capillaryPressureModel(const dictionary& dict)
    :
        dict_(dict)
    {}

    virtual ~capillaryPressureModel() = default;

    // ---- Interface ----
    virtual void correct
    (
        volScalarField& pcow,
        volScalarField& pcgo,
        volScalarField& dpcow_dsw,
        volScalarField& dpcow_dsg,
        volScalarField& dpcgo_dsw,
        volScalarField& dpcgo_dsg,
        const volScalarField& Sa,
        const volScalarField& Sb,
        const volScalarField& zeroField
    ) const = 0;

    // Factory
    static autoPtr<capillaryPressureModel> New
    (
        const dictionary& transportProperties
    );
};

} // End namespace Foam

#endif
