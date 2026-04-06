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
    Foam::capillaryPressureModel::ashooriPC

Description
    Ashoori capillary pressure model for two-phase flow.

    Reference:
    \verbatim
        Ashoori, E., Marchesin, D., & Rossen, W. R. (2010, September).
        Roles of transient and local equilibrium foam behavior in porous media–traveling wave.
        ECMOR XII-12th European Conference on the Mathematics of Oil Recovery (pp. cp-163). European Association of Geoscientists & Engineers.
    \endverbatim

    The default model coefficients are
    \verbatim
        ashooriPC
        {
            Swc   0.2;
            Sgr   0.18;
            pc0   330;
            c     0.01;
        }
    \endverbatim

SourceFiles
    ashooriPC.C

\*---------------------------------------------------------------------------*/

#ifndef ashooriPC_H
#define ashooriPC_H

#include "capillaryPressureModel.H"

namespace Foam
{

class ashooriPC
:
    public capillaryPressureModel
{
    scalar Swc_;
    scalar Sgr_;
    scalar pc0_;
    scalar c_;

public:
    TypeName("ashooriPC");

    ashooriPC(const dictionary& dict);

    virtual void correct
    (
        volScalarField& pc,
        volScalarField& dpcds,
        const volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
