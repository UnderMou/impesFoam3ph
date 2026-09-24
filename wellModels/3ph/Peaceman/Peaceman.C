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

List<well> Peaceman::readWells(const dictionary& wellsDict)
{
    wordList wellNames = wellsDict.toc();

    List<well> wells(wellNames.size());

    rateFns_.setSize(wellNames.size());
    bhpFns_.setSize(wellNames.size());
    FaInjFns_.setSize(wellNames.size());
    FbInjFns_.setSize(wellNames.size());
    CsInjFns_.setSize(wellNames.size());

    forAll(wellNames, i)
    {
        const word& name = wellNames[i];

        const dictionary& wellDict = wellsDict.subDict(name);

        wells[i].name = name;

        wells[i].cells = labelList(wellDict.lookup("cells"));

        wells[i].injector = readBool(wellDict.lookup("injector"));

        const bool rateGiven = wellDict.found("rate");
        const bool bhpGiven = wellDict.found("bhp");

        if (rateGiven == bhpGiven)
        {
            FatalIOErrorInFunction(wellDict)
                << "Well " << name << " must give exactly one of rate or bhp"
                << exit(FatalIOError);
        }

        wells[i].bhpControl = bhpGiven;

        wells[i].rate   = 0.0;
        wells[i].bhp    = 0.0;
        wells[i].Fa_inj = 0.0;
        wells[i].Fb_inj = 0.0;
        wells[i].Cs_inj = 0.0;

        // Evaluated by updateWellInputs() at each time step
        if (wells[i].bhpControl)
        {
            bhpFns_.set(i, Function1<scalar>::New("bhp", wellDict));
        }
        else
        {
            rateFns_.set(i, Function1<scalar>::New("rate", wellDict));
        }

        if (wells[i].injector)
        {
            FaInjFns_.set(i, Function1<scalar>::New("Fa_inj", wellDict));
            FbInjFns_.set(i, Function1<scalar>::New("Fb_inj", wellDict));
            CsInjFns_.set(i, Function1<scalar>::New("Cs_inj", wellDict));
        }

        wells[i].rhoWell = 0.0;

        wells[i].radius = wellDict.lookupOrDefault<scalar>("radius", 0.497); 

        Info << "Reading well: " << name << nl;
    }

    return wells;
}

Peaceman::Peaceman(const dictionary& dict)
:
    wellModel(dict),
    vertDir_(dict.lookupOrDefault<vector>("vertDir", vector(0, 0, 1)))
{
    wells_ = readWells(wellsDict_);
    Info << "Vertical direction: " << vertDir_ << nl << endl;
}

void Peaceman::source_pEqn
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
) 
{
    wellCoeff = scalar(0.0);
    wellSource = scalar(0.0);
    const scalarField& V = p.mesh().V();
    const volVectorField& C = p.mesh().C();

    // Correct rhoWell value
    forAll(wells_, w)
    {
        well& well = wells_[w];

        // producer well
        if (!well.injector)
        {
            scalar QaProd = 0.0;
            scalar QbProd = 0.0;
            scalar QcProd = 0.0;

            forAll(well.cells, j)
            {
                const label celli = well.cells[j];

                const scalar Qt = qt[celli]*V[celli];
                const scalar Qa = qa[celli]*V[celli];
                const scalar Qb = qb[celli]*V[celli];
                const scalar Qc = Qt - Qa - Qb;
 
                QaProd += -min(Qa, scalar(0.0));
                QbProd += -min(Qb, scalar(0.0));
                QcProd += -min(Qc, scalar(0.0));
            }

            const scalar QProd = QaProd + QbProd + QcProd;

            if (QProd > SMALL)
            {
                well.rhoWell = (rho_a*QaProd + rho_b*QbProd + rho_c*QcProd)/QProd;
            }
            else
            {
                // No previous production:
                // use following initialization value

                scalar numerator   = 0.0;
                scalar denominator = 0.0;

                forAll(well.cells, j)
                {
                    const label celli = well.cells[j];

                    numerator += WI[celli]*(mob_a[celli]*rho_a + mob_b[celli]*rho_b + mob_c[celli]*rho_c);
                    denominator += WI[celli]*mob_t[celli];
                }

                well.rhoWell = numerator/denominator;
            }

            Info<< "Well " << well.name
                << ": QaProd = " << QaProd
                << ", QbProd = " << QbProd
                << ", QcProd = " << QcProd
                << ", rhoWell = " << well.rhoWell
                << nl;
        }
        // injection well
        else
        {
            well.rhoWell = well.Fa_inj*rho_a + well.Fb_inj*rho_b + (1.0 - well.Fa_inj - well.Fb_inj)*rho_c;
        }
    }

    // Evaluate well contributions to pressure equation
    forAll(wells_, w)
    {
        const well& well = wells_[w];
        const scalar ratePerPerf = well.rate/scalar(well.cells.size()); // Qi = Q/Nperf

        // Info << "name:" << well.name << ", injector?" << well.injector << ", bhp controlled?" << well.bhpControl << endl;
        // Info << "control values: bhp=" << well.bhp << ", rate=" << well.rate << endl;
        // Info << well.name << ": C[" << well.cells[0] << "] = " << C[well.cells[0]] << endl;

        forAll(well.cells, j)
        {
            label celli = well.cells[j];

            // Info << celli << " " << endl;
            if (wells_[w].bhpControl)
            {
                // q = WI*mob*(p_bh - p - rhoWell*g*Delta_z) --- prescribed bhp
                wellCoeff[celli] += WI[celli]*mob_t[celli]/V[celli];
                wellSource[celli] += WI[celli]*mob_t[celli]*(wells_[w].bhp - well.rhoWell * (g.value() & (C[well.cells[0]] - C[celli])))/V[celli];
            }
            else
            {
                // q = Q/Nperf --- prescribed rate
                wellSource[celli] += ratePerPerf/V[celli]; 
            }
        }
    }
    pEqn += fvm::Sp(wellCoeff, p);
    pEqn -= wellSource;
}

void Peaceman::source_SaEqn
(
    fvScalarMatrix& SaEqn,
    const volScalarField& Sa, 
    const volScalarField& Fa, 
    const volScalarField& p,
    scalar t,
    volScalarField& qa
)
{   
    Info<< "qa: min = " << gMin(qa.internalField()) << ", max = " << gMax(qa.internalField()) << nl << endl;
    SaEqn -= qa;
}

void Peaceman::source_SbEqn
(
    fvScalarMatrix& SbEqn,
    const volScalarField& Sb, 
    const volScalarField& Fb, 
    const volScalarField& p,
    scalar t,
    volScalarField& qb
)
{   
    Info<< "qb: min = " << gMin(qb.internalField()) << ", max = " << gMax(qb.internalField()) << nl << endl;
    SbEqn -= qb;
}

void Peaceman::correct
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
) 
{
    const volVectorField& C = p.mesh().C();
    const scalarField& V = qt.mesh().V();

    // Correct qt 
    qt = scalar(0.0);

    forAll(wells_, w)
    {
        const well& well = wells_[w];

        if (well.bhpControl)
        {
            forAll(well.cells, j)
            {
                label celli = well.cells[j];
                
                // q = WI*mob*(p_bh - p - rhoWell*g*Delta_z) --- prescribed bhp
                qt[celli] += WI[celli]*mob_t[celli]*(well.bhp - p[celli] - well.rhoWell * (g.value() & (C[well.cells[0]] - C[celli])))/V[celli]; 
            }
        }
        else
        {
            const scalar ratePerPerf = well.rate/scalar(well.cells.size());

            forAll(well.cells, j)
            {
                label celli = well.cells[j];
                
                // q = Q/Nperf --- prescribed rate
                qt[celli] += ratePerPerf/V[celli];                
            }
        }
    
    }

    // check total rate: Sum(qt*V) == 0
    scalar totalRate = 0.0;
    forAll(qt.internalField(), celli)
    {
        totalRate += qt[celli]*V[celli];
    }
    Info<< "Integrated total well rate = " << totalRate << nl << endl;

    // Correct bhp and rate
    forAll(wells_, w)
    {
        well& well = wells_[w];

        if (!well.bhpControl)
        {
            scalar sumT = 0.0;
            scalar sumTp = 0.0;

            forAll(well.cells, j)
            {
                const label celli = well.cells[j];

                const scalar Ti = WI[celli]*mob_t[celli];

                sumT  += Ti;
                sumTp += Ti*p[celli] + Ti*well.rhoWell * (g.value() & (C[well.cells[0]] - C[celli]));
            }

            well.bhp = (well.rate + sumTp)/sumT;

            forAll(well.cells, j)
            {
                const label celli = well.cells[j];

                // well pressure = BHP - rhoWell*g*Delta_z
                p_bh[celli] = well.bhp - well.rhoWell * (g.value() & (C[well.cells[0]] - C[celli]));
            }
        }
        else
        {
            well.rate = 0.0;

            forAll(well.cells, j)
            {
                const label celli = well.cells[j];

                well.rate += WI[celli]*mob_t[celli]*(well.bhp - p[celli] - well.rhoWell * (g.value() & (C[well.cells[0]] - C[celli])));    
                
                // well pressure = BHP - rhoWell*g*Delta_z
                p_bh[celli] = well.bhp - well.rhoWell * (g.value() & (C[well.cells[0]] - C[celli]));
            }
            Info << "Well " << well.name << " bhp = " << well.bhp << ", rate = " << well.rate << nl << endl;
        }
        
    }

    // Correct qa
    qa = scalar(0.0);

    forAll(wells_, w)
    {
        const well& well = wells_[w];

        forAll(well.cells, j)
        {
            label celli = well.cells[j];

            if (well.injector)
            {
                qa[celli] += well.Fa_inj*qt[celli];              
            }
            else
            {
                qa[celli] += Fa[celli]*qt[celli];
            }
        }
    }
        
    // Correct qb
    qb = scalar(0.0);

    forAll(wells_, w)
    {
        const well& well = wells_[w];

        forAll(well.cells, j)
        {
            label celli = well.cells[j];

            if (well.injector)
            {
                qb[celli] += well.Fb_inj*qt[celli];              
            }
            else
            {
                qb[celli] += Fb[celli]*qt[celli];
            }
        }
    }

    // Correct qs
    qs = scalar(0.0);

    forAll(wells_, w)
    {
        const well& well = wells_[w];
        
        if (well.injector)
        {
            forAll(well.cells, j)
            {
                label celli = well.cells[j];
                qs[celli] += well.Cs_inj*qb[celli]; 
            }
        }
        else
        {
            forAll(well.cells, j)
            {
                label celli = well.cells[j];
                qs[celli] += Cs[celli]*qb[celli]; 
            }
        }
    }

}

void Peaceman::initialize
(
    volScalarField& WI,
    const volScalarField& K
)
{
    const scalarField& V = WI.mesh().V();

    scalarField re = Foam::sqrt(V/M_PI);

    Info<< "min(re) = " << gMin(re) << " max(re) = " << gMax(re) << nl << endl;

    forAll(wells_, w)
    {
        const well& well = wells_[w];

        forAll(well.cells, j)
            {
                label celli = well.cells[j];

                // Evaluate cell height
                const pointField& points = WI.mesh().points();
                const cell& c = WI.mesh().cells()[celli];
                const labelList pointLabels = c.labels(WI.mesh().faces());
                scalar sMin = GREAT;
                scalar sMax = -GREAT;

                forAll(pointLabels, pi)
                {
                    const point& pt = points[pointLabels[pi]];

                    const scalar s = pt & vertDir_;

                    sMin = min(sMin, s);
                    sMax = max(sMax, s);
                }

                const scalar cellHeight = sMax - sMin;
                Info<< "Well " << well.name << ": cellHeight = " << cellHeight << nl << endl;

                // Evaluate Well Index
                if (re[celli]/well.radius <= 1.0)
                {
                    FatalErrorInFunction
                        << "Invalid Peaceman geometry: re/rw <= 1 for well "
                        << well.name << ", cell " << celli
                        << exit(FatalError);
                }
                WI[celli] = 2.0*M_PI*K[celli]*cellHeight / Foam::log(re[celli]/well.radius); // TODO: implement skin factor correction
            }
        Info<< "Well " << well.name << ": radius = " << well.radius << nl << endl;
    }

    Info<< "min(WI) = " << gMin(WI.internalField())
        << " max(WI) = " << gMax(WI.internalField())
        << nl << endl;

    Info<< "WI OK" << nl << endl;

} 

} // End namespace Foam