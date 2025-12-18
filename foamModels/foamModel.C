#include "foamModel.H"

namespace Foam
{

defineTypeNameAndDebug
(
    foamModel,
    0
);

defineRunTimeSelectionTable
(
    foamModel,
    dictionary
);

autoPtr<foamModel>
foamModel::New
(
    const dictionary& transportProperties,
    foamAuxFields* aux
)
{
    word modelName
    (
        transportProperties.lookup("foamModel")
    );

    Info<< "Selecting foam model: "
        << modelName << nl << endl;

    const dictionary& modelDict =
        transportProperties.subDict(modelName);

    auto iter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (iter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown foam model "
            << modelName << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<foamModel>
    (
        iter()(modelDict, aux)
    );
}


} // End namespace Foam
