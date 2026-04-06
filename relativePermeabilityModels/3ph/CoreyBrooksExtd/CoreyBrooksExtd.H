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
    Foam::relativePermeabilityModel::CoreyBrooksExtd

Description
    Corey-Brooks relative permeability model for three-phase flow.

    Obs.: The difference with Foam::relativePermeabilityModel::CoreyBrooks
    is that this model allows for phase 'i' saturation bellow residuals (Si < Sri) by setting kri = 0
    and Si > 1 - Srj for j != i by setting kri = kri_max.
    
    Reference:
    \verbatim
        Lyu, X., Voskov, D., Tang, J., & Rossen, W. R. (2021). 
        Simulation of foam enhanced-oil-recovery processes using operator-based linearization approach.
        SPE Journal, 26(04), 2287-2304.
    \endverbatim

    The default model coefficients are
    \verbatim
        CoreyBrooksExtd
        {
            Sa_min    0.0;
            Sb_min    0.0;
            Sc_min    0.0;
            a_exp     2.0;
            b_exp     2.0;
            c_exp     2.0;
            kra_max   1.0;
            krb_max   1.0;
            krc_max   1.0;
        }
    \endverbatim

SourceFiles
    CoreyBrooksExtd.C

\*---------------------------------------------------------------------------*/

#ifndef CoreyBrooksExtd_H
#define CoreyBrooksExtd_H

#include "relativePermeabilityModel.H"

namespace Foam
{

class CoreyBrooksExtd
:
    public relativePermeabilityModel
{
    scalar Sa_min_;
    scalar Sb_min_;
    scalar Sc_min_;
    scalar a_exp_;
    scalar b_exp_;
    scalar c_exp_;
    scalar kra_max_;
    scalar krb_max_;
    scalar krc_max_;

public:
    TypeName("CoreyBrooksExtd");

    CoreyBrooksExtd(const dictionary& dict);

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
