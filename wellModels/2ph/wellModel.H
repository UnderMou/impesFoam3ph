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
    Foam::wellModel

Description
    This class supplies a well model template, 
    and is the base class for a number of well models implementation.

SourceFiles
    wellModelModel.C

\*---------------------------------------------------------------------------*/

#ifndef wellModel_H
#define wellModel_H

#include "fvCFD.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "Function1.H"
#include "wellStruct.H"

namespace Foam
{

class wellModel
{
protected:
    const dictionary& dict_;
    const dictionary& wellsDict_;
    List<well> wells_;

    // Time-dependent inputs of each well
    PtrList<Function1<scalar>> rateFns_;
    PtrList<Function1<scalar>> bhpFns_;
    PtrList<Function1<scalar>> FbInjFns_;
    PtrList<Function1<scalar>> CsInjFns_;

    virtual List<well> readWells(const dictionary& wellsDict) = 0;

    // Check if all wells_ are rate controlled, 
    // then check for sum(rates) = 0 condition
    // (incompressible closed system).
    void checkRateBalance()
    {
        bool allRateControlled = true;

        forAll(wells_, w)
        {
            if (wells_[w].bhpControl)
            {
                allRateControlled = false;
                break;
            }
        }

        if (allRateControlled)
        {
            scalar Qtotal = 0.0;

            forAll(wells_, w)
            {
                Qtotal += wells_[w].rate;
            }

            if (mag(Qtotal) > SMALL)
            {
                FatalErrorInFunction
                    << "All wells are rate controlled, but Sum(Q) != 0." << nl
                    << "Sum(Q) = " << Qtotal << nl
                    << exit(FatalError);
            }
        }
    }

public:
    TypeName("wellModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        wellModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    wellModel(const dictionary& dict)
    :
        dict_(dict),
        wellsDict_(dict.subDict("wells")),
        wells_()
    {
    }

    virtual ~wellModel() = default;

    // Evaluate the well inputs at time t
    virtual void updateWellInputs(const scalar t)
    {
        forAll(wells_, w)
        {
            well& wl = wells_[w];

            Info<< "Well " << wl.name << ": ";

            if (wl.bhpControl)
            {
                wl.bhp = bhpFns_[w].value(t);
                Info<< "bhp = " << wl.bhp;
            }
            else
            {
                wl.rate = rateFns_[w].value(t);
                Info<< "rate = " << wl.rate;
            }

            if (wl.injector)
            {
                wl.Fb_inj = FbInjFns_[w].value(t);
                wl.Cs_inj = CsInjFns_[w].value(t);
                Info<< ", Fb_inj = " << wl.Fb_inj << ", Cs_inj = " << wl.Cs_inj;

                if (wl.Fb_inj < -SMALL || wl.Fb_inj > 1 + SMALL)
                {
                    FatalErrorInFunction
                        << "Well " << wl.name << ": Fb_inj = " << wl.Fb_inj
                        << " at t = " << t << " is outside [0, 1]"
                        << exit(FatalError);
                }
            }

            Info<< endl;
        }

        checkRateBalance();
    }

    // ---- Interface ----

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
        const dimensionedVector& g,
        const volScalarField& qt,
        const volScalarField& qb
    ) = 0;

    virtual void source_SbEqn
    (
        fvScalarMatrix& SbEqn,
        const volScalarField& Sb,
        const volScalarField& Fb,
        const volScalarField& p,
        scalar t,
        volScalarField& qb
    ) = 0;

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
    ) = 0;

    virtual void initialize
    (
        volScalarField& WI,
        const volScalarField& K
    )
    {
        // Default: do nothing
    }
    
    // Factory
    static autoPtr<wellModel> New
    (
        const dictionary& wellProperties
    );
};

} // End namespace Foam

#endif
