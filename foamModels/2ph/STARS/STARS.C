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

#include <cmath>
#include "STARS.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    STARS,
    0
);

addToRunTimeSelectionTable
(
    foamModel,
    STARS,
    dictionary
);

STARS::STARS(const dictionary& dict, foamAuxFields* aux)
:
    foamModel(dict, aux),
    fmmob_(readScalar(dict.lookup("fmmob"))),
    SF_(readScalar(dict.lookup("SF"))),
    sfbet_(readScalar(dict.lookup("sfbet"))),
    fmcap_(readScalar(dict.lookup("fmcap"))),
    epcap_(readScalar(dict.lookup("epcap"))),
    sigma_ba_(readScalar(dict.lookup("sigma_ba"))),
    fmsurf_(readScalar(dict.lookup("fmsurf"))),
    epsurf_(readScalar(dict.lookup("epsurf"))),
    FshearActive_(dict.parent().lookupOrDefault<Switch>("FshearActive",false)),
    FsurfActive_(dict.parent().lookupOrDefault<Switch>("FsurfActive",false)),
    ncaModel_(dict.parent().lookupOrDefault<word>("ncaModel","zeng"))
{}

void STARS::correct_Fdry
(
    volScalarField& Fdry,
    const volScalarField Sb
) const
{
    Fdry = 0.5 + (1.0 / M_PI) * Foam::atan(sfbet_ * (Sb - SF_));
}

void STARS::correct_Nca
(
    volScalarField& Nca,
    const volVectorField U,
    const volScalarField& K,
    const volVectorField& gradP
) const
{
    if (ncaModel_ == "zeng") {
        Nca = (mu_b_.value() * mag(U)) / sigma_ba_;
    } else if (ncaModel_ == "boeije"){
        Nca = (K * mag(gradP)) / sigma_ba_; 
    } else {
        FatalErrorInFunction
            << "Unknown Nca model "
            << ncaModel_ << nl
            << exit(FatalError);
    }
    Info << "Nca model: " << ncaModel_ << endl;
}

void STARS::correct_Fshear
(
    volScalarField& Fshear,
    volScalarField& Nca
) const
{
    // Correct at cell center
    forAll(Nca,i){
        if (Nca[i] >= fmcap_){
            Fshear[i] = Foam::pow(fmcap_/Nca[i], epcap_);
        } else {
            Fshear[i] = 1.0;
        }
    } 

    // Correct at the boundaries
    forAll(Nca.boundaryField(), patchi)
    {
        fvPatchScalarField& NcaPatch = Nca.boundaryFieldRef()[patchi];
        fvPatchScalarField& FshearPatch = Fshear.boundaryFieldRef()[patchi];

        forAll(NcaPatch, facei)
        {
            if (NcaPatch[facei] >= fmcap_){
                FshearPatch[facei] = Foam::pow(fmcap_/NcaPatch[facei], epcap_);
            } else {
                FshearPatch[facei] = 1.0;
            }
        }
    }
}

void STARS::correct_Fsurf
(
    volScalarField& Fsurf,
    volScalarField& Cs
) const
{
    // Correct at cell center
    forAll(Cs,i){
        if (Cs[i] < fmsurf_){
            Fsurf[i] = Foam::pow(Cs[i]/fmsurf_, epsurf_);
        } else {
            Fsurf[i] = 1.0;
        }
    } 

    // Correct at the boundaries
    forAll(Cs.boundaryField(), patchi)
    {
        fvPatchScalarField& CsPatch = Cs.boundaryFieldRef()[patchi];
        fvPatchScalarField& FsurfPatch = Fsurf.boundaryFieldRef()[patchi];

        forAll(CsPatch, facei)
        {
            if (CsPatch[facei] < fmsurf_){
                FsurfPatch[facei] = Foam::pow(CsPatch[facei]/fmsurf_, epsurf_);
            } else {
                FsurfPatch[facei] = 1.0;
            }
        }
    }
}

void STARS::correct_MRF
(
    volScalarField& MRF,
    const volScalarField Fdry,
    const volScalarField Fshear,
    const volScalarField Fsurf
) const
{
    MRF = 1.0 / (1.0 + fmmob_*Fdry*Fshear*Fsurf);
}

void STARS::correct
(
    volScalarField& kra,
    const volVectorField U,
    const volScalarField& Sa,
    const volScalarField& Sb,
    const surfaceScalarField& phia,
    const volScalarField& eps,
    const volScalarField& K,
    const volVectorField& gradP
) const
{

    if (!aux_ || !aux_->Fdry)
    {
        FatalErrorInFunction
            << "STARS requires auxiliary fields but none provided"
            << exit(FatalError);
    }

    volScalarField& Nca    = *aux_->Nca;
    volScalarField& Fdry   = *aux_->Fdry;
    volScalarField& Fshear = *aux_->Fshear;
    volScalarField& Fsurf  = *aux_->Fsurf;
    volScalarField& MRF    = *aux_->MRF;
    volScalarField& Cs     = *aux_->Cs;

    correct_Fdry(Fdry, Sb);
    if (FshearActive_)
    {
        correct_Nca(Nca, U, K, gradP);
        correct_Fshear(Fshear, Nca);
    }
    if (FsurfActive_)
    {
        correct_Fsurf(Fsurf, Cs);
    }
    correct_MRF(MRF, Fdry, Fshear, Fsurf);
    kra *= MRF;
}

} // End namespace Foam
