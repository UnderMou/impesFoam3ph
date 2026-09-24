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
    Constant injection rate model for three-phase flow.

    Reference:
    \verbatim

    \endverbatim

    The default model coefficients are
    \verbatim

        wellModel ConstantInjection;

        ConstantInjection
        {
            vertDir (0 1 0);
            
            wells
            {
                injector1
                {
                name        INJ1;
                cells       (0 1 2 3 4 5 6 7 8 9);  // cell index
                injector    true;                   // false = producer
                rate        1e-4;                   // injection rate [m3/s]
                Fa_inj      0.0;                    // injection phase a action
                Fb_inj      0.0;                    // injection phase b action
                Cs_inj      0.0;                    // injection phase concentration [kg/m3]
                }

                producer1
                {
                name        PROD1;
                cells       (90 91 92 93 94 95 96 97 98 99); 
                injector    false;
                rate        -1e-4;                  // negative for a producer
                }
            }
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
        const scalar& rho_c,
        const volScalarField& mob_a,
        const volScalarField& mob_b,
        const volScalarField& mob_c,
        const dimensionedVector& g,
        const volScalarField& qt,
        const volScalarField& qa,
        const volScalarField& qb
    ) override;

    virtual void source_SaEqn
    (
        fvScalarMatrix& SaEqn,
        const volScalarField& Sa,
        const volScalarField& Fa,
        const volScalarField& p,
        scalar t,
        volScalarField& qa
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
        volScalarField& qa,
        volScalarField& qb,
        const volScalarField& Fa,
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
        const scalar& rho_c,
        const volScalarField& mob_a,
        const volScalarField& mob_b,
        const volScalarField& mob_c,
        const dimensionedVector& g
    ) override;
};

} // End namespace Foam

#endif
