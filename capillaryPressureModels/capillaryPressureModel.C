#include "capillaryPressureModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    capillaryPressureModel,
    0
);

defineRunTimeSelectionTable
(
    capillaryPressureModel,
    dictionary
);

autoPtr<capillaryPressureModel>
capillaryPressureModel::New
(
    const dictionary& transportProperties
)
{
    word modelName
    (
        transportProperties.lookup("capillaryPressureModel")
    );

    Info<< "Selecting capillary pressure model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        transportProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown capillary pressure model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<capillaryPressureModel>
    (
        iter()(modelDict)
    );
}


} // End namespace Foam
