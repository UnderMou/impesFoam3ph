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
    Foam::foamModel::ashoori

Description
    Ashoori foam model for two-phase flow.

    Reference:
    \verbatim
        Ashoori, E., Marchesin, D., & Rossen, W. R. (2010, September).
        Roles of transient and local equilibrium foam behavior in porous media–traveling wave.
        ECMOR XII-12th European Conference on the Mathematics of Oil Recovery (pp. cp-163). European Association of Geoscientists & Engineers.
    \endverbatim

    The default model coefficients are
    \verbatim

        ashoori
        {
        Kc      0.01;
        A       400;
        Cmrf    18500;
        Sw_ast  0.37;
        nMAX    8.0e13;
        }

    \endverbatim

SourceFiles
    ashoori.C

\*---------------------------------------------------------------------------*/

#ifndef ashoori_H
#define ashoori_H

#include "fvCFD.H"
#include "foamModel.H"

namespace Foam
{

class ashoori
:
    public foamModel
{
    scalar Cmrf_;
    scalar Kc_;
    scalar nMAX_;
    scalar Sw_ast_;
    scalar A_;

public:
    TypeName("ashoori");

    ashoori(const dictionary& dict, foamAuxFields* aux);

    // ashoori functions
    void correct_nD_LE
    (
        volScalarField& nDLE,
        volScalarField Sb
    ) const;

    void correct_rgrc
    (
        volScalarField& rgrc,
        volScalarField& nDLE,
        const volScalarField Sb,
        const volScalarField nD
    ) const;

    void correct_MRF
    (
        volScalarField& MRF,
        const volScalarField nD
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
