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

    forAll(wellNames, i)
    {
        const word& name = wellNames[i];

        const dictionary& wellDict =
            wellsDict.subDict(name);

        wells[i].name = name;

        wells[i].cells =
            labelList(wellDict.lookup("cells"));

        wells[i].injector =
            readBool(wellDict.lookup("injector"));

        wells[i].bhpControl =
            readBool(wellDict.lookup("bhpControl"));

        wells[i].target =
            readScalar(wellDict.lookup("target"));
        
        if (wells[i].bhpControl)
        {
            wells[i].bhp  = wells[i].target;
            wells[i].rate = 0.0;
        }
        else
        {
            wells[i].bhp  = 0.0;
            wells[i].rate = wells[i].target;
        }
        
        if (wells[i].injector)
        {
            wells[i].Fb_inj = readScalar(wellDict.lookup("Fb_inj"));;
            wells[i].Cs_inj = readScalar(wellDict.lookup("Cs_inj"));
        }
        else
        {
            wells[i].Fb_inj = 0.0;
            wells[i].Cs_inj = 0.0;
        }

        Info << "Reading well: " << name << nl;
    }

    return wells;
}

Peaceman::Peaceman(const dictionary& dict)
:
    wellModel(dict)
{
    wells_ = readWells(wellsDict_);
    checkRateBalance();
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
    const volScalarField& mob_a,
    const volScalarField& mob_b,
    const dimensionedVector& g
) 
{
    wellCoeff = scalar(0.0);
    wellSource = scalar(0.0);

    const scalarField& V = p.mesh().V();

    forAll(wells_, w)
    {
        const well& well = wells_[w];
        const scalar ratePerPerf = well.rate/scalar(well.cells.size());

        Info << "name:" << well.name << ", injector?" << well.injector << ", bhp controlled?" << well.bhpControl << endl;
        Info << "control values: bhp=" << well.bhp << ", rate=" << well.rate << endl;

        forAll(well.cells, j)
        {
            label celli = well.cells[j];

            // Info << celli << " " << endl;
            
            if (wells_[w].bhpControl)
            {
                // q = WI*mob*(p_bh - p)/V
                wellCoeff[celli] += WI[celli]*mob_t[celli]/V[celli];
                wellSource[celli] += WI[celli]*mob_t[celli]*wells_[w].bhp/V[celli];
            }
            else
            {
                // q = prescribed rate
                wellSource[celli] += ratePerPerf/V[celli]; // wells_[w].rate;
            }
        }
    }
    pEqn += fvm::Sp(wellCoeff, p);
    pEqn -= wellSource;
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
) 
{
    // Correct qt and p_bh
    qt = scalar(0.0);

    const scalarField& V = qt.mesh().V();

    forAll(wells_, w)
    {
        const well& well = wells_[w];

        if (well.bhpControl)
        {
            forAll(well.cells, j)
            {
                label celli = well.cells[j];
                
                // q = WI*mob*(p_bh - p)/V
                qt[celli] += WI[celli]*mob_t[celli]*(well.bhp - p[celli])/V[celli];  
            }
        }
        else
        {
            const scalar ratePerPerf = well.rate/scalar(well.cells.size());

            forAll(well.cells, j)
            {
                label celli = well.cells[j];
                
                // q = Q/V
                qt[celli] += ratePerPerf/V[celli];                
            }
        }
    
    }

    // check total rate
    scalar totalRate = 0.0;

    forAll(qt.internalField(), celli)
    {
        totalRate += qt[celli]*V[celli];
    }

    Info<< "Integrated total well rate = "
        << totalRate << nl << endl;

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
                sumTp += Ti*p[celli];
            }

            well.bhp = (well.rate + sumTp)/sumT;

            forAll(well.cells, j)
            {
                const label celli = well.cells[j];

                p_bh[celli] = well.bhp;
            }
        }
        else
        {
            well.rate = 0.0;

            forAll(well.cells, j)
            {
                const label celli = well.cells[j];

                well.rate += WI[celli]*mob_t[celli]*(well.bhp - p[celli]);    
            
                p_bh[celli] = well.bhp;
            }
            Info << "Well " << well.name << " bhp = " << well.bhp << ", rate = " << well.rate << nl << endl;
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

} // End namespace Foam
