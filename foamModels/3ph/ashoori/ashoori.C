/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "ashoori.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    ashoori,
    0
);

addToRunTimeSelectionTable
(
    foamModel,
    ashoori,
    dictionary
);

ashoori::ashoori(const dictionary& dict, foamAuxFields* aux)
:
    foamModel(dict, aux),
    Cmrf_(readScalar(dict.lookup("Cmrf"))),
    Kc_(readScalar(dict.lookup("Kc"))),
    nMAX_(readScalar(dict.lookup("nMAX"))),
    Sw_ast_(readScalar(dict.lookup("Sw_ast"))),
    A_(readScalar(dict.lookup("A")))
{}

void ashoori::correct_nD_LE
(
    volScalarField& nDLE,
    volScalarField Sb
) const
{
    // Correct at cell center
    forAll(nDLE,i){
        if (Sb[i] > Sw_ast_){
            nDLE[i] = Foam::tanh(A_ * (Sb[i] - Sw_ast_));
        }
    }

    // Correct at the boundaries
    forAll(Sb.boundaryField(), patchi)
    {
        fvPatchScalarField& nDLEPatch = nDLE.boundaryFieldRef()[patchi];
        fvPatchScalarField& SbPatch = Sb.boundaryFieldRef()[patchi];

        forAll(SbPatch, facei)
        {
            if (SbPatch[facei] > Sw_ast_){
                nDLEPatch[facei] = Foam::tanh(A_ * (SbPatch[facei] - Sw_ast_));
            }
        }
    }
}

void ashoori::correct_rgrc
(
    volScalarField& rgrc,
    volScalarField& nDLE,
    const volScalarField Sb,
    const volScalarField nD
) const
{
    correct_nD_LE(nDLE, Sb);
    rgrc = Kc_ * nMAX_ * (nDLE - nD);
}

void ashoori::correct_MRF
(
    volScalarField& MRF,
    const volScalarField nD
) const
{
    MRF = Cmrf_ * nD + scalar(1.0);
}

void ashoori::correct
(
    volScalarField& kra,
    const volVectorField U,
    const volScalarField& Sa,
    const volScalarField& Sb,
    volScalarField& Sc,
    const surfaceScalarField& phia,
    const volScalarField& eps,
    const volScalarField& K,
    const volVectorField& gradP
) const
{

    if (!aux_ || !aux_->rgrc)
    {
        FatalErrorInFunction
            << "ashoori requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& nD        = *aux_->nD;
    volScalarField& rgrc      = *aux_->rgrc;
    volScalarField& F         = *aux_->F;
    volScalarField& AcumCoeff = *aux_->AcumCoeff_nD;
    volScalarField& nDLE      = *aux_->nDLE;
    volScalarField& MRF       = *aux_->MRF;

    F = -fvc::div(phia); // + qa // TODO
    AcumCoeff = eps*(Sa + VSMALL);
    correct_rgrc(rgrc, nDLE, Sb, nD);
    
    // nD equation
    fvScalarMatrix nDEqn
    (
        AcumCoeff*fvm::ddt(nD) + fvc::div(phia, nD) 
        ==
        (eps*(Sa + VSMALL)/nMAX_)*rgrc + fvc::Sp(-F, nD) // + qf // TODO
    );
    nDEqn.solve();

    correct_MRF(MRF, nD);
    
    kra /= MRF;
}

} // End namespace Foam
