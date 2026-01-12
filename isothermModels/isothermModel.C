#include "isothermModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    isothermModel,
    0
);

defineRunTimeSelectionTable
(
    isothermModel,
    dictionary
);

autoPtr<isothermModel>
isothermModel::New
(
    const dictionary& transportProperties,
    foamAuxFields* aux
)
{
    word modelName
    (
        transportProperties.lookup("isothermModel")
    );

    Info<< "Selecting isotherm model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        transportProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown isotherm model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<isothermModel>
    (
        iter()(modelDict, aux)
    );
}


} // End namespace Foam
