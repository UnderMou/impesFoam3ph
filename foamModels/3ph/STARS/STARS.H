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
    Foam::foamModel::STARS

Description
    STARS foam model for three-phase flow.

    Reference:
    \verbatim
        Lyu, X., Voskov, D., Tang, J., & Rossen, W. R. (2021). 
        Simulation of foam enhanced-oil-recovery processes using operator-based linearization approach.
        SPE Journal, 26(04), 2287-2304.
    \endverbatim

    The default model coefficients are
    \verbatim

        foamModel STARS;
        FsurfActive false;
        FshearActive false;
        FoilActive true;

        STARS
        {
        fmmob     2000.0;
        SF        0.3;
        sfbet     32000.0;
        fmcap     1.0;
        epcap     0.0;
        sigma_ba  1.0;
        fmsurf    1.0;
        epsurf    0.0;
        fmoil     0.3;
        floil     0.1;
        epoil     3.0;
        }

    \endverbatim

SourceFiles
    STARS.C

\*---------------------------------------------------------------------------*/

#ifndef STARS_H
#define STARS_H

#include "foamModel.H"

namespace Foam
{

class STARS
:
    public foamModel
{
    scalar fmmob_;
    scalar SF_;
    scalar sfbet_;
    scalar fmcap_;
    scalar epcap_;
    scalar sigma_ba_;
    scalar fmsurf_;
    scalar epsurf_;
    scalar fmoil_;
    scalar floil_;
    scalar epoil_;

    Switch FshearActive_;
    Switch FsurfActive_;
    Switch FoilActive_;

    word   ncaModel_;

public:
    TypeName("STARS");

    STARS(const dictionary& dict, foamAuxFields* aux);

    // STARS Fi functions
    void correct_Fdry
    (
        volScalarField& Fdry,
        const volScalarField& Sb
    ) const;

    void correct_Nca
    (
        volScalarField& Nca,
        const volVectorField& U,
        const volScalarField& K,
        const volVectorField& gradP
    ) const;

    void correct_Fshear
    (
        volScalarField& Fshear,
        volScalarField& Nca
    ) const;

    void correct_Fsurf
    (
        volScalarField& Fsurf,
        volScalarField& Cs
    ) const;

    void correct_Foil
    (
        volScalarField& Foil,
        volScalarField& Sc
    ) const;

    void correct_MRF
    (
        volScalarField& MRF,
        const volScalarField& Fdry,
        const volScalarField& Fshear,
        const volScalarField& Fsurf,
        const volScalarField& Foil
    ) const;

    // Correct
    virtual void correct
    (
        volScalarField& kra,
        const volVectorField U,
        const volScalarField& Sa,
        const volScalarField& Sb,
        volScalarField& Sc,
        const surfaceScalarField& phia,
        const volScalarField& eps,
        const volScalarField& K,
        const volVectorField& gradP
    ) const;
};

} // End namespace Foam

#endif
