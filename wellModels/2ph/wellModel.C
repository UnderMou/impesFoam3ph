#include "wellModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    wellModel,
    0
);

defineRunTimeSelectionTable
(
    wellModel,
    dictionary
);

autoPtr<wellModel>
wellModel::New
(
    const dictionary& wellProperties
)
{
    word modelName
    (
        wellProperties.lookup("wellModel")
    );

    Info<< "Selecting well model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        wellProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown well model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<wellModel>
    (
        iter()(modelDict)
    );
}


} // End namespace Foam
