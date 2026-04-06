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
    Sc_min_(readScalar(dict.lookup("Sc_min"))),
    a_exp_(readScalar(dict.lookup("a_exp"))),
    b_exp_(readScalar(dict.lookup("b_exp"))),
    c_exp_(readScalar(dict.lookup("c_exp"))),
    kra_max_(readScalar(dict.lookup("kra_max"))),
    krb_max_(readScalar(dict.lookup("krb_max"))),
    krc_max_(readScalar(dict.lookup("krc_max")))
{}

void CoreyBrooksExtd::correct
(
    volScalarField& kra,
    volScalarField& krb,
    volScalarField& krc,
    volScalarField& Sa,
    volScalarField& Sb,
    volScalarField& Csw
) const
{

    // Evaluate oil saturation
    volScalarField Sc = scalar(1.0) - Sa - Sb;

    // Effective saturations (normalization)
    volScalarField Se_a = (Sa - Sa_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_);
    volScalarField Se_b = (Sb - Sb_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_);
    volScalarField Se_c = (Sc - Sc_min_) / (1.0 - Sa_min_ - Sb_min_ - Sc_min_);

    // Correct kra
    // cell centered values
    forAll(kra,i){
        if (Sa[i]<=Sa_min_ + SMALL){
            kra[i]=0.0;
        }
        else if (Sa[i]>=1.0 - Sb_min_ - Sc_min_ - SMALL){
            kra[i]=kra_max_;
        }
        else{
            kra[i] = kra_max_ * Foam::pow(Se_a[i], a_exp_);
        }
    }
    // face centered values
    forAll(kra.boundaryField(), patchi)
    {
        fvPatchScalarField& kraPatch = kra.boundaryFieldRef()[patchi];
        fvPatchScalarField& SaPatch = Sa.boundaryFieldRef()[patchi];
        fvPatchScalarField& Se_aPatch = Se_a.boundaryFieldRef()[patchi];

        forAll(kraPatch, facei)
        {
            if (SaPatch[facei]<=Sa_min_ + SMALL){
                kraPatch[facei]=0.0;
            }
            else if (SaPatch[facei]>=1.0 - Sb_min_ - Sc_min_ - SMALL){
                kraPatch[facei]=kra_max_;
            }
            else{
                kraPatch[facei] = kra_max_ * Foam::pow(Se_aPatch[facei], a_exp_);
            }
        }
    }

    // Correct krb
    // cell centered values
    forAll(krb,i){
        if (Sb[i]<=Sb_min_ + SMALL){
            krb[i]=0.0;
        }
        else if (Sb[i]>=1.0 - Sa_min_ - Sc_min_ - SMALL){
            krb[i]=krb_max_;
        }
        else{
            krb[i] = krb_max_ * Foam::pow(Se_b[i], b_exp_);
        }
    }
    // face centered values
    forAll(krb.boundaryField(), patchi)
    {
        fvPatchScalarField& krbPatch = krb.boundaryFieldRef()[patchi];
        fvPatchScalarField& SbPatch = Sb.boundaryFieldRef()[patchi];
        fvPatchScalarField& Se_bPatch = Se_b.boundaryFieldRef()[patchi];

        forAll(krbPatch, facei)
        {
            if (SbPatch[facei]<=Sb_min_ + SMALL){
                krbPatch[facei]=0.0;
            }
            else if (SbPatch[facei]>=1.0 - Sa_min_ - Sc_min_ - SMALL){
                krbPatch[facei]=krb_max_;
            }
            else{
                krbPatch[facei] = krb_max_ * Foam::pow(Se_bPatch[facei], b_exp_);
            }
        }
    }

    // Correct krc
    // cell centered values
    forAll(krc,i){
        if (Sc[i]<=Sc_min_ + SMALL){
            krc[i]=0.0;
        }
        else if (Sc[i]>=1.0 - Sb_min_ - Sa_min_ - SMALL){
            krc[i]=krc_max_;
        }
        else{
            krc[i] = krc_max_ * Foam::pow(Se_c[i], c_exp_);
        }
    }
    // face centered values
    forAll(krc.boundaryField(), patchi)
    {
        fvPatchScalarField& krcPatch = krc.boundaryFieldRef()[patchi];
        fvPatchScalarField& ScPatch = Sc.boundaryFieldRef()[patchi];
        fvPatchScalarField& Se_cPatch = Se_c.boundaryFieldRef()[patchi];

        forAll(krcPatch, facei)
        {
            if (ScPatch[facei]<=Sc_min_ + SMALL){
                krcPatch[facei]=0.0;
            }
            else if (ScPatch[facei]>=1.0 - Sb_min_ - Sa_min_ - SMALL){
                krcPatch[facei]=krc_max_;
            }
            else{
                krcPatch[facei] = krc_max_ * Foam::pow(Se_cPatch[facei], c_exp_);
            }
        }
    }
}

} // End namespace Foam
