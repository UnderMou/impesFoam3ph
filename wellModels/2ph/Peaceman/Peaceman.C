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

#include "Peaceman.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug
(
    Peaceman,
    0
);

addToRunTimeSelectionTable
(
    wellModel,
    Peaceman,
    dictionary
);

Peaceman::Peaceman(const dictionary& dict)
:
    wellModel(dict),
    Fb_inj_(readScalar(dict.lookup("Fb_inj"))),
    inj_maxBHP_(readScalar(dict.lookup("inj_maxBHP"))),
    inj_maxRate_(readScalar(dict.lookup("inj_maxRate"))),
    prod_maxBHP_(readScalar(dict.lookup("prod_maxBHP"))),
    prod_maxRate_(readScalar(dict.lookup("prod_maxRate"))),
    slugTime_w_(readScalar(dict.lookup("slugTime_w"))),
    slugTime_g_(readScalar(dict.lookup("slugTime_g"))),
    CsInj_(readScalar(dict.lookup("CsInj"))),
    t0_(dict.lookupOrDefault<scalar>("t0", 0.0))
    
{
    cycleTime_ = slugTime_w_ + slugTime_g_;
    inj_bhp_control_ = false;
    prod_bhp_control_ = false;
    Info << "t0_ = " << t0_ << endl;
}

tmp<fvScalarMatrix> Peaceman::implicitSource_pEqn
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
) const
{
    if (inj_bhp_control_)
    {
        // BHP is specified, calculate injection rate
        return fvm::Sp((WI * mob_t)*isWell, p);
    }
    else
    {
        // Injection rate is specified, calculate BHP
        return tmp<fvScalarMatrix>(new fvScalarMatrix(p, dimless));
    }
}

tmp<volScalarField> Peaceman::explicitSource_pEqn
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
) const
{   
    if (inj_bhp_control_)
    {
        // BHP is specified, calculate injection rate
        return (WI * mob_t * p_bh) * isWell;
    }
    else
    {
        // Injection rate is specified, calculate BHP
        return qt_inj - qt_prod;
    }
}

// virtual tmp<fvScalarMatrix> implicitSource_SEqn
// (
//     volScalarField& qb,
//     const volScalarField& Sb, 
//     const volScalarField& Fb, 
//     const volScalarField& p,
//     const volScalarField& qt_inj,
//     const volScalarField& qt_prod,
//     scalar t
// ) const
// {
//     if (bhp_control_)
//     {
//         // BHP is specified, calculate injection rate
//         // TODO:
//     }
//     else
//     {
//         // Injection rate is specified, calculate BHP
//         // TODO:
//     }
// } 

tmp<volScalarField> Peaceman::explicitSource_SEqn
(
    volScalarField& qb,
    const volScalarField& Sb, 
    const volScalarField& Fb, 
    const volScalarField& p,
    const volScalarField& qt_inj,
    const volScalarField& qt_prod,
    scalar t
) const
{
    
    scalar t_cycle = std::fmod(t - 1e7, cycleTime_);

    if (t_cycle < 0.0)
    {
        t_cycle += cycleTime_;
    }

    bool inWater = (t_cycle >= 0.0 && t_cycle < slugTime_w_);
    bool inGas   = (t_cycle >= slugTime_w_ && t_cycle < cycleTime_);

    Info << "Time: " << t << " - In water slug: " << inWater << " - In gas slug: " << inGas << endl;

    if(inWater)
    {
        qb = Fb_inj_ * qt_inj - Fb * qt_prod; 
    }
    
    if(inGas)
    {
        qb = (scalar(1.0) - Fb_inj_) * qt_inj - Fb * qt_prod; 
    }

    return qb;
}

void Peaceman::correct
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
) 
{
    if (inj_bhp_control_)
    {
        // BHP is specified, calculate injection rate
        qt = (WI * mob_t * (p_bh - p)) * isWell;

        // scalar qtMax = gMax(qt.internalField());
        // qt = qtMax * isWell;

        qt_inj = Foam::mag(max(qt, scalar(0)));
        qt_prod = Foam::mag(min(qt, scalar(0)));

        // evaluate BHP
        p_bh = (p + (qt / (WI * mob_t))) * isWell;
        p_bh = min(p_bh, inj_maxBHP_); // limit BHP to maximum allowed value
    }
    else
    {
        // Injection rate is specified, calculate BHP
        p_bh = (p + (qt / (WI * mob_t))) * isWell;
    }

    // check if the BHP exceeds the maximum allowed value
    if (gMax(p_bh.internalField()) >= inj_maxBHP_)
    {
        Info << "BHP = " << gMax(p_bh.internalField()) << ", exceeds maximum allowed value. Switching to BHP control." << endl;
        inj_bhp_control_ = true;

        p_bh = min(p_bh, inj_maxBHP_); // limit BHP to maximum allowed value
    }
    
    Info << "BHP = " << gMax(p_bh) << endl;
    Info << "well constraint: " << (inj_bhp_control_ ? "BHP control" : "Rate control") << endl;

    // Correct qs
    qs = CsInj_*qb*pos(qb) + Cs*qb*neg(qb);
}

} // End namespace Foam
