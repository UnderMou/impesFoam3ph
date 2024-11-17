/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "codedFvModelTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
#line 34 "/home/anderson/OpenFOAM/anderson-9/run/impes_case_v3/constant/fvModels/SwSource"

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 4d18e215335a88a92567a544bee3526e1935bd16
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void SwSource_4d18e215335a88a92567a544bee3526e1935bd16(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SwSourceFvModelscalarSource, 0);

addRemovableToRunTimeSelectionTable
(
    fvModel,
    SwSourceFvModelscalarSource,
    dictionary
);


const char* const SwSourceFvModelscalarSource::SHA1sum =
    "4d18e215335a88a92567a544bee3526e1935bd16";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SwSourceFvModelscalarSource::
SwSourceFvModelscalarSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(coeffs(), mesh)
{
    if (false)
    {
        Info<<"construct SwSource sha1: 4d18e215335a88a92567a544bee3526e1935bd16"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SwSourceFvModelscalarSource::
~SwSourceFvModelscalarSource()
{
    if (false)
    {
        Info<<"destroy SwSource sha1: 4d18e215335a88a92567a544bee3526e1935bd16\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SwSourceFvModelscalarSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"SwSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 38 "/home/anderson/OpenFOAM/anderson-9/run/impes_case_v3/constant/fvModels/SwSource"
// Pout<< "**codeAddSup**" << endl;
        // const Time& time = mesh().time();
        // const scalarField& V = mesh().V();
        // scalarField& heSource = eqn.source();
        // heSource -= 0.1*V; //*sqr(time.value())*V;

        const Time& time = mesh().time();
        const scalar currentTime = time.value();

        const scalarField& V = mesh().V();
        scalarField& SwSource = eqn.source();
        
        // // =========== BOX1 SELECTION ========================
        // // ===================================================
        // // Retrieve the cell set named "box1"
        // cellSet selectedCells(mesh(), "box1");
        // // Iterate over the cell set and apply the source term
        // forAllConstIter(labelHashSet, selectedCells, iter)
        // {
        //     label cellI = iter.key();
        //     // SwSource[cellI] -= (0.1) * V[cellI];
            
        //     const scalar x = mesh().C()[cellI].component(0);
        //     const scalar y = 0.5;

        //     SwSource[cellI] -= ((x-pow(x,2))) * V[cellI];
        // }



        // =========== ALL SELECTION ========================
        // ===================================================
        forAll(SwSource, i)
        {
            const scalar x = mesh().C()[i].component(0);
            const scalar y = 0.5;

            // SwSource[i] -= ((x-pow(x,2))*(y-pow(y,2)) + 
            //                 2.0*currentTime*(x-pow(x,2)+y-pow(y,2))) * V[i];

            SwSource[i] -= (x/4.0 - pow(x,2)/4.0 + 0.5*currentTime) * V[i];
        };
//}}} end code
}


void SwSourceFvModelscalarSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"SwSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 86 "/home/anderson/OpenFOAM/anderson-9/run/impes_case_v3/constant/fvModels/SwSource"
Pout<< "**codeAddRhoSup**" << endl;
//}}} end code
}


void SwSourceFvModelscalarSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"SwSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 91 "/home/anderson/OpenFOAM/anderson-9/run/impes_case_v3/constant/fvModels/SwSource"
Pout<< "**codeAddAlphaRhoSup**" << endl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace fv

// ************************************************************************* //

