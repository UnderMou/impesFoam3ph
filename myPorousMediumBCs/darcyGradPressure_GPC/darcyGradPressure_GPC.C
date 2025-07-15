/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

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

#include "darcyGradPressure_GPC.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyGradPressure_GPC::darcyGradPressure_GPC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_("Mf"),
    phiName_("phi"),
    phiGfName_("phiG"),
    phiPcName_("phiPc")
{}

Foam::darcyGradPressure_GPC::darcyGradPressure_GPC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    phiGfName_(dict.lookupOrDefault<word>("phiG","phiG")),
    phiPcName_(dict.lookupOrDefault<word>("phiPc","phiPc"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::darcyGradPressure_GPC::darcyGradPressure_GPC
(
    const darcyGradPressure_GPC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_),
    phiPcName_(ptf.phiPcName_)
{}

// Foam::darcyGradPressure_GPC::darcyGradPressure_GPC
// (
//     const darcyGradPressure_GPC& ptf
// )
//     :
//     fixedGradientFvPatchScalarField(ptf),
//     MfName_(ptf.MfName_),
//     phiName_(ptf.phiName_)
// {}

Foam::darcyGradPressure_GPC::darcyGradPressure_GPC
(
    const darcyGradPressure_GPC& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_),
    phiPcName_(ptf.phiPcName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyGradPressure_GPC::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& Mf=
        patch().lookupPatchField<surfaceScalarField, scalar>(MfName_);

    const fvsPatchField<scalar>& phi=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvsPatchField<scalar>& phiGf=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiGfName_);
    
    const fvsPatchField<scalar>& phiPc=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiPcName_);

    gradient() = - (phi-phiGf-phiPc)/(Mf+ROOTVSMALL)/(patch().magSf());

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void Foam::darcyGradPressure_GPC::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "phiG", "phiG", phiGfName_);
    writeEntryIfDifferent<word>(os, "phiPc", "phiPc", phiPcName_);
    writeEntry(os, "Mf", MfName_);
    writeEntry(os, "phi", phiName_);
    writeEntry(os, "phiG", phiGfName_);
    writeEntry(os, "phiPc", phiPcName_);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyGradPressure_GPC
);
}

// ************************************************************************* //
