#include "relativePermeabilityModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    relativePermeabilityModel,
    0
);

defineRunTimeSelectionTable
(
    relativePermeabilityModel,
    dictionary
);

autoPtr<relativePermeabilityModel>
relativePermeabilityModel::New
(
    const dictionary& transportProperties
)
{
    word modelName
    (
        transportProperties.lookup("relativePermeabilityModel")
    );

    Info<< "Selecting relative permeability model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        transportProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown relative permeability model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<relativePermeabilityModel>
    (
        iter()(modelDict)
    );
}


} // End namespace Foam
