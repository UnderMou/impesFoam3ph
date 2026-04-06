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
    Foam::relativePermeabilityModel::vanGenuchten

Description
    vanGenuchten relative permeability model for two-phase flow.

    Reference:
    \verbatim
        Horgue, P., Soulaine, C., Franc, J., Guibert, R., & Debenest, G. (2015).
        An open-source toolbox for multiphase flow in porous media.
        Computer Physics Communications, 187, 217-226.
    \endverbatim

    The default model coefficients are
    \verbatim
        vanGenuchten
        {
            Sa_min    0.0;
            Sb_min    0.0;
            n         2.0;
            m         0.5;
            kra_max   1.0;
            krb_max   1.0;
        }
    \endverbatim

SourceFiles
    vanGenuchten.C

\*---------------------------------------------------------------------------*/

#ifndef vanGenuchten_H
#define vanGenuchten_H

#include "relativePermeabilityModel.H"

namespace Foam
{

class vanGenuchten
:
    public relativePermeabilityModel
{
    scalar Sa_min_;
    scalar Sb_min_;
    scalar n_;
    scalar m_;
    scalar kra_max_;
    scalar krb_max_;

public:
    TypeName("vanGenuchten");

    vanGenuchten(const dictionary& dict);

    virtual void correct
    (
        volScalarField& kra,
        volScalarField& krb,
        volScalarField& Sb,
        volScalarField& Csw
    ) const;
};

} // End namespace Foam

#endif
