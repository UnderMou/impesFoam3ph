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
    Foam::wellModel::ConstantInjection

Description
    Constant injection rate model for two-phase flow.

    Reference:
    \verbatim
        Horgue, P., Soulaine, C., Franc, J., Guibert, R., & Debenest, G. (2015).
        An open-source toolbox for multiphase flow in porous media.
        Computer Physics Communications, 187, 217-226.
    \endverbatim

    The default model coefficients are
    \verbatim

        wellModel ConstantInjection;

        ConstantInjection
        {
            // co-injection
            Fb_inj    0.3;
        }

    \endverbatim

SourceFiles
    ConstantInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ConstantInjection_H
#define ConstantInjection_H

#include "wellModel.H"

namespace Foam
{

class ConstantInjection
:
    public wellModel
{


protected:
    virtual List<well> readWells
    (
        const dictionary& wellsDict
    )override;

public:
    TypeName("ConstantInjection");

    ConstantInjection(const dictionary& dict);

    virtual void source_pEqn
    (
        fvScalarMatrix& pEqn,
        const volScalarField& p,
        const volScalarField& mob_t,
        const volScalarField& WI,
        volScalarField& wellCoeff,
        volScalarField& wellSource,
        const scalar& rho_a,
        const scalar& rho_b,
        const volScalarField& mob_a,
        const volScalarField& mob_b,
        const dimensionedVector& g
        
    ) override;

    virtual void source_SbEqn
    (
        fvScalarMatrix& SbEqn,
        const volScalarField& Sb, 
        const volScalarField& Fb, 
        const volScalarField& p,
        scalar t,
        volScalarField& qb
    ) override;

    virtual void correct
    (
        volScalarField& qt,
        volScalarField& qb,
        const volScalarField& Fb,
        const volScalarField& p, 
        scalar t,
        const volScalarField& mob_t,
        const volScalarField& WI,
        volScalarField& p_bh,
        volScalarField& qs,
        const volScalarField& Cs,
        const scalar& rho_a,
        const scalar& rho_b,
        const volScalarField& mob_a,
        const volScalarField& mob_b,
        const dimensionedVector& g
    ) override;
};

} // End namespace Foam

#endif
