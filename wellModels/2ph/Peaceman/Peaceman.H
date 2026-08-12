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
    Foam::wellModel::Peaceman

Description
    Peaceman well model operates according to the following constraints:
        > Bottom hole pressure (BHP) is specified -> injection rate is calculated
        > Injection rate is specified -> BHP is calculated

    The default model coefficients are
    \verbatim

        wellModel Peaceman;

        Peaceman
        {
            maxBHP 1000000; // [Pa] Maximum BHP for injection well
            maxRate 1e-3;   // [m3/s /m3 = 1/s] Maximum injection rate            
        }

    \endverbatim

SourceFiles
    Peaceman.C

\*---------------------------------------------------------------------------*/

#ifndef Peaceman_H
#define Peaceman_H

#include "fvCFD.H"
#include "wellModel.H"

namespace Foam
{

class Peaceman
:
    public wellModel
{
    scalar Fb_inj_;

    // peaceman
    scalar inj_maxBHP_;
    scalar inj_maxRate_;
    bool inj_bhp_control_;
    scalar prod_maxBHP_;
    scalar prod_maxRate_;
    bool prod_bhp_control_;

    // slug
    scalar slugTime_w_;
    scalar slugTime_g_;
    scalar cycleTime_;
    scalar t_cycle_;
    bool inWater_; 
    bool inGas_; 

    // surfactant
    scalar CsInj_;

    scalar t0_;

public:
    TypeName("Peaceman");

    Peaceman(const dictionary& dict);

    virtual tmp<fvScalarMatrix> implicitSource_pEqn
    (
        volScalarField& qb,
        const volScalarField& Sb, 
        const volScalarField& Fb, 
        const volScalarField& p,
        const volScalarField& qt_inj,
        const volScalarField& qt_prod,
        scalar t,
        const volScalarField& mob_t,
        const volScalarField& WI,
        const volScalarField& p_bh,
        const volScalarField& isWell
    ) const override;

    virtual tmp<volScalarField> explicitSource_pEqn
    (
        volScalarField& qb,
        const volScalarField& Sb, 
        const volScalarField& Fb, 
        const volScalarField& p,
        const volScalarField& qt_inj,
        const volScalarField& qt_prod,
        scalar t,
        const volScalarField& mob_t,
        const volScalarField& WI,
        const volScalarField& p_bh,
        const volScalarField& isWell
    ) const override;

    // virtual tmp<fvScalarMatrix> implicitSource_SEqn
    // (
    //     volScalarField& qb,
    //     const volScalarField& Sb, 
    //     const volScalarField& Fb, 
    //     const volScalarField& p,
    //     const volScalarField& qt_inj,
    //     const volScalarField& qt_prod,
    //     scalar t
    // ) const override;

    virtual tmp<volScalarField> explicitSource_SEqn
    (
        volScalarField& qb,
        const volScalarField& Sb, 
        const volScalarField& Fb, 
        const volScalarField& p,
        const volScalarField& qt_inj,
        const volScalarField& qt_prod,
        scalar t
    ) const override;

    virtual void correct
    (
        volScalarField& qt,
        volScalarField& qb,
        const volScalarField& Fb,
        const volScalarField& p, 
        scalar t,
        volScalarField& qt_inj,
        volScalarField& qt_prod,
        const volScalarField& mob_t,
        const volScalarField& WI,
        volScalarField& p_bh,
        const volScalarField& isWell,
        volScalarField& qs,
        const volScalarField& Cs
    ) override;

};

} // End namespace Foam

#endif
