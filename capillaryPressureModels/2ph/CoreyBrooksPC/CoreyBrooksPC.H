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
    Foam::capillaryPressureModel::CoreyBrooksPC

Description
    Corey-Brooks capillary pressure model for two-phase flow.

    Reference:
    \verbatim
        Horgue, P., Soulaine, C., Franc, J., Guibert, R., & Debenest, G. (2015).
        An open-source toolbox for multiphase flow in porous media.
        Computer Physics Communications, 187, 217-226.
    \endverbatim

    The default model coefficients are
    \verbatim
        CoreyBrooksPC
        {
            Spc_irr 0.0;
            Spc_max 1.0;
            pc0     1000.0;
            alpha   0.5;
        }
    \endverbatim

SourceFiles
    CoreyBrooksPC.C

\*---------------------------------------------------------------------------*/

#ifndef CoreyBrooksPC_H
#define CoreyBrooksPC_H

#include "capillaryPressureModel.H"

namespace Foam
{

class CoreyBrooksPC
:
    public capillaryPressureModel
{
    scalar Spc_irr_;
    scalar Spc_max_;
    scalar pc0_;
    double alpha_;

public:
    TypeName("CoreyBrooksPC");

    CoreyBrooksPC(const dictionary& dict);

    virtual void correct
    (
        volScalarField& pc,
        volScalarField& dpcds,
        const volScalarField& Sb
    ) const;
};

} // End namespace Foam

#endif
