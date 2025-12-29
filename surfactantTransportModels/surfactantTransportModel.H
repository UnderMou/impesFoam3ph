#ifndef surfactantTransportModel_H
#define surfactantTransportModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "foamAuxFields.H"

namespace Foam
{

class surfactantTransportModel
{
protected:
    const dictionary& dict_;
    foamAuxFields* aux_;

    Switch active_;
    dimensionedScalar rho_b_;

public:
    TypeName("surfactantTransportModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        surfactantTransportModel,
        dictionary,
        (
            const dictionary& dict,
            foamAuxFields* aux
        ),
        (dict, aux)
    );

    surfactantTransportModel(const dictionary& dict, foamAuxFields* aux)
    :
        dict_(dict),
        aux_(aux),
        active_(dict.parent().lookupOrDefault<Switch>("surfTranspActive",false)),
        rho_b_(dict.parent().lookup("rho_b"))
    {
        if (!active_)
        {
            Info<< "Surfactant concentration transport disabled\n" << endl;
            return;
        }

        Info<< "Surfactant concentration transport enabled\n" << endl;
    }

    virtual ~surfactantTransportModel() = default;

    // ---- Interface ----
    Switch active() const { return active_; }

    virtual void correct
    (
        const volScalarField& Sb,
        const surfaceScalarField& phib,
        const volScalarField& eps
    ) const = 0;

    // Factory
    static autoPtr<surfactantTransportModel> New
    (
        const dictionary& transportProperties,
        foamAuxFields* aux
    );
};

} // End namespace Foam

#endif
