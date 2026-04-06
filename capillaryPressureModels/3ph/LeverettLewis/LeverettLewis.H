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

Class
    Foam::capillaryPressureModel::LeverettLewis

Description
    LeverettLewis capillary pressure model for three-phase flow.

    Reference:
    \verbatim
        Abreu, E., Furtado, F., & Pereira, F. (2004).
        On the Numerical Simulation of Three‐Phase Reservoir Transport Problems. 
        Transport Theory and Statistical Physics, 33(5-7), 503-526.
    \endverbatim

    The default model coefficients are
    \verbatim
        LeverettLewis
        {
            epsilon     0.005;
        }
    \endverbatim

SourceFiles
    LeverettLewis.C

\*---------------------------------------------------------------------------*/


#ifndef LeverettLewis_H
#define LeverettLewis_H

#include "capillaryPressureModel.H"

namespace Foam
{

class LeverettLewis
:
    public capillaryPressureModel
{
    scalar epsilon_;

public:
    TypeName("LeverettLewis");

    LeverettLewis(const dictionary& dict);

    virtual void correct
    (
        volScalarField& pcow,
        volScalarField& pcgo,
        volScalarField& dpcow_dsw,
        volScalarField& dpcow_dsg,
        volScalarField& dpcgo_dsw,
        volScalarField& dpcgo_dsg,
        const volScalarField& Sa,
        const volScalarField& Sb,
        const volScalarField& zeroField
    ) const;
};

} // End namespace Foam

#endif
