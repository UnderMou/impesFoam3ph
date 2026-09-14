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
#include "wellStruct.H"

namespace Foam
{

class wellModel
{
protected:
    const dictionary& dict_;
    const dictionary& wellsDict_;
    List<well> wells_;

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
        const scalar& rho_c,
        const volScalarField& mob_a,
        const volScalarField& mob_b,
        const volScalarField& mob_c,
        const dimensionedVector& g,
        const volScalarField& qt,
        const volScalarField& qa,
        const volScalarField& qb
    ) = 0;

    virtual void source_SaEqn
    (
        fvScalarMatrix& SaEqn,
        const volScalarField& Sa,
        const volScalarField& Fa,
        const volScalarField& p,
        scalar t,
        volScalarField& qa
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
