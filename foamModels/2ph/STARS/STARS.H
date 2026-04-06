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
    STARS foam model for two-phase flow.

    Reference:
    \verbatim
        Valdez, A. R., Rocha, B. M., da Fonseca Façanha, J. M., de Souza, A. V. O., Pérez-Gramatges, A., Chapiro, G., & Santos, R. W. D. (2022). 
        Foam-assisted water–gas flow parameters: from core-flood experiment to uncertainty quantification and sensitivity analysis. 
        Transport in Porous Media, 144(1), 189-209.
    \endverbatim

    The default model coefficients are
    \verbatim

        foamModel STARS;
        FshearActive false;
        FsurfActive false;

        STARS
        {
        fmmob     269.342055;
        SF        0.43766393;
        sfbet     787.457944;
        fmcap     1.0;
        epcap     0.0;
        sigma_ba  1.0; 
        fmsurf    1.0;
        epsurf    0.0;
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

    Switch FshearActive_;
    Switch FsurfActive_;

    word   ncaModel_;

public:
    TypeName("STARS");

    STARS(const dictionary& dict, foamAuxFields* aux);

    // STARS Fi functions
    void correct_Fdry
    (
        volScalarField& Fdry,
        const volScalarField Sb
    ) const;

    void correct_Nca
    (
        volScalarField& Nca,
        const volVectorField U,
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

    void correct_MRF
    (
        volScalarField& MRF,
        const volScalarField Fdry,
        const volScalarField Fshear,
        const volScalarField Fsurf
    ) const;

    // Correct
    virtual void correct
    (
        volScalarField& kra,
        const volVectorField U,
        const volScalarField& Sa,
        const volScalarField& Sb,
        const surfaceScalarField& phia,
        const volScalarField& eps,
        const volScalarField& K,
        const volVectorField& gradP
    ) const;
};

} // End namespace Foam

#endif
