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

#include "fvMatrices.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class wellModel
{
protected:
    const dictionary& dict_;

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
        dict_(dict)
    {}

    virtual ~wellModel() = default;

    // ---- Interface ----
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
    ) const 
    {
        return tmp<fvScalarMatrix>(new fvScalarMatrix(p, dimless));
    }

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
    ) const 
    {
        tmp<volScalarField> tSource(new volScalarField(p));
        tSource.ref() = dimensionedScalar("zero", dimless, 0.0);
        return tSource;
    }

    virtual tmp<fvScalarMatrix> implicitSource_SEqn
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
        return tmp<fvScalarMatrix>(new fvScalarMatrix(Sb, dimless));
    }

    virtual tmp<volScalarField> explicitSource_SEqn
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
        tmp<volScalarField> tSource(new volScalarField(Sb));
        tSource.ref() = dimensionedScalar("zero", dimless, 0.0);
        return tSource;
    }

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
    ) 
    {
        // Default implementation does nothing,
        // derived classes can override this method 
        // to apply corrections to qt based on the well model.
    }
    
    // Factory
    static autoPtr<wellModel> New
    (
        const dictionary& wellProperties
    );
};

} // End namespace Foam

#endif
