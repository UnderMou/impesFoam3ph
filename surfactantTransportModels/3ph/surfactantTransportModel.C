#include "surfactantTransportModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    surfactantTransportModel,
    0
);

defineRunTimeSelectionTable
(
    surfactantTransportModel,
    dictionary
);

autoPtr<surfactantTransportModel>
surfactantTransportModel::New
(
    const dictionary& transportProperties,
    foamAuxFields* aux
)
{
    word modelName
    (
        transportProperties.lookup("surfactantTransportModel")
    );

    Info<< "Selecting surfactant transport model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        transportProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown surfactant transport model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<surfactantTransportModel>
    (
        iter()(modelDict, aux)
    );
}


} // End namespace Foam
