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

#include "ashooriPC.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    ashooriPC,
    0
);

addToRunTimeSelectionTable
(
    capillaryPressureModel,
    ashooriPC,
    dictionary
);

ashooriPC::ashooriPC(const dictionary& dict)
:
    capillaryPressureModel(dict),
    Swc_(readScalar(dict.lookup("Swc"))),
    Sgr_(readScalar(dict.lookup("Sgr"))),
    pc0_(readScalar(dict.lookup("pc0"))),
    c_(readScalar(dict.lookup("c")))
{

}

void ashooriPC::correct
(
    volScalarField& pc,
    volScalarField& dpcds,
    const volScalarField& Sb
) const
{   
    // correct pc and its derivative
    dpcds = ((pc0_ * c_ * (-scalar(1.0) + Swc_ + Sgr_)) / Foam::pow(Sb - Swc_ + VSMALL, 2))*
                            Foam::pow((scalar(1.0) - Sb - Sgr_)/(Sb - Swc_),c_-1); 
                            
    pc = pc0_ * Foam::pow(scalar(1.0) - Sb - Sgr_, c_) / (Sb - Swc_ + VSMALL); 
}

} // End namespace Foam
