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
    Foam::relativePermeabilityModel::Mehrabi

Description
    Mehrabi relative permeability model for three-phase flow.
    It accounts for oil mobilization due to surfactant, i.e., Sor = Sor(Cs)
    and krw_max = krw_max(Cs).

    Reference:
    \verbatim
        Mehrabi, M., Sepehrnoori, K., & Delshad, M. (2022).
        Displacement theory of low-tension gas flooding.
        Transport in Porous Media, 142(3), 475-491.
    \endverbatim

    The default model coefficients are
    \verbatim
        Mehrabi
        {
            Sa_min    0.25;
            Sb_min    0.1;
            Sc_min_1  0.0;
            Sc_min_2  0.3;
            a_exp     2.0;
            b_exp     2.0;
            c_exp     2.0;
            kra_max   1.0;
            krb_max_1 0.8;
            krb_max_2 0.3;
            krc_max   1.0;

            Cr        1e-2;
        }
    \endverbatim

SourceFiles
    Mehrabi.C

\*---------------------------------------------------------------------------*/

#ifndef Mehrabi_H
#define Mehrabi_H

#include "relativePermeabilityModel.H"

namespace Foam
{

class Mehrabi
:
    public relativePermeabilityModel
{
    scalar Sa_min_;
    scalar Sb_min_;
    scalar Sc_min_1_;
    scalar Sc_min_2_;
    scalar a_exp_;
    scalar b_exp_;
    scalar c_exp_;
    scalar kra_max_;
    scalar krb_max_1_;
    scalar krb_max_2_;
    scalar krc_max_;
    scalar Cr_;

public:
    TypeName("Mehrabi");

    Mehrabi(const dictionary& dict);

    scalar correctSor
    (
        scalar Csw
    ) const;

    scalar correctKrw0
    (
        scalar Csw
    ) const;

    virtual void correct
    (
        volScalarField& kra,
        volScalarField& krb,
        volScalarField& krc,
        volScalarField& Sa,
        volScalarField& Sb,
        volScalarField& Csw
    ) const;
};

} // End namespace Foam

#endif
