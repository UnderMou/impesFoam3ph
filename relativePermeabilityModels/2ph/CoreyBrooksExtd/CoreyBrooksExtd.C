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

#include "CoreyBrooksExtd.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    CoreyBrooksExtd,
    0
);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    CoreyBrooksExtd,
    dictionary
);

CoreyBrooksExtd::CoreyBrooksExtd(const dictionary& dict)
:
    relativePermeabilityModel(dict),
    Sa_min_(readScalar(dict.lookup("Sa_min"))),
    Sb_min_(readScalar(dict.lookup("Sb_min"))),
    a_exp_(readScalar(dict.lookup("a_exp"))),
    b_exp_(readScalar(dict.lookup("b_exp"))),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_(readScalar(dict.lookup("krb_max")))
{}

void CoreyBrooksExtd::correct
(
    volScalarField& kra,
    volScalarField& krb,
    volScalarField& Sb,
    volScalarField& Csw
) const
{

    // Effective saturations (normalization)
    volScalarField Se = (Sb - Sb_min_) / (1.0 - Sa_min_ - Sb_min_);
    volScalarField Sa = scalar(1.0) - Sb;

    // Correct kra
    // cell centered values
    forAll(kra,i){
        if (Sa[i]<=Sa_min_ + SMALL){
            kra[i]=0.0;
        }
        else if (Sa[i]>=1.0 - Sb_min_ - SMALL){
            kra[i]=kra_max_;
        }
        else{
            kra[i] = kra_max_ * Foam::pow(scalar(1.0) - Se[i], a_exp_);
        }
    }
    // face centered values
    forAll(kra.boundaryField(), patchi)
    {
        fvPatchScalarField& kraPatch = kra.boundaryFieldRef()[patchi];
        fvPatchScalarField& SaPatch = Sa.boundaryFieldRef()[patchi];
        fvPatchScalarField& SePatch = Se.boundaryFieldRef()[patchi];

        forAll(kraPatch, facei)
        {
            if (SaPatch[facei]<=Sa_min_ + SMALL){
                kraPatch[facei]=0.0;
            }
            else if (SaPatch[facei]>=1.0 - Sb_min_ - SMALL){
                kraPatch[facei]=kra_max_;
            }
            else{
                kraPatch[facei] = kra_max_ * Foam::pow(scalar(1.0) - SePatch[facei], a_exp_);
            }
        }
    }

    // Correct krb
    // cell centered values
    forAll(krb,i){
        if (Sb[i]<=Sb_min_ + SMALL){
            krb[i]=0.0;
        }
        else if (Sb[i]>=1.0 - Sa_min_ - SMALL){
            krb[i]=krb_max_;
        }
        else{
            krb[i] = krb_max_ * Foam::pow(Se[i], b_exp_);
        }
    }
    // face centered values
    forAll(krb.boundaryField(), patchi)
    {
        fvPatchScalarField& krbPatch = krb.boundaryFieldRef()[patchi];
        fvPatchScalarField& SbPatch = Sb.boundaryFieldRef()[patchi];
        fvPatchScalarField& SePatch = Se.boundaryFieldRef()[patchi];

        forAll(krbPatch, facei)
        {
            if (SbPatch[facei]<=Sb_min_ + SMALL){
                krbPatch[facei]=0.0;
            }
            else if (SbPatch[facei]>=1.0 - Sa_min_ - SMALL){
                krbPatch[facei]=krb_max_;
            }
            else{
                krbPatch[facei] = krb_max_ * Foam::pow(SePatch[facei], b_exp_);
            }
        }
    }
}

} // End namespace Foam
