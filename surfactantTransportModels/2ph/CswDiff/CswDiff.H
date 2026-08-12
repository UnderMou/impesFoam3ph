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
    Foam::surfactantTransportModel::CswDiff

Description
    Surfactant transport by aqueous phase (only) model for two-phase flow accounting with surfactant
    concentration diffusion.

    Reference:
    \verbatim
        CHUNG, F. TH. Modeling of surfactant transport and adsorption in porous media. National Inst. for Petroleum and Energy Research, Bartlesville, OK (USA), 1991.
    \endverbatim

    The default model coefficients are
    \verbatim

        CswDiff
        {
            D 1.0e-9;
            rho_s 2.0e3;

            isothermModel Henry;

            Henry
            {
                Kd    0.05;
            }

        }

    \endverbatim

SourceFiles
    CswDiff.C

\*---------------------------------------------------------------------------*/

#ifndef CswDiff_H
#define CswDiff_H

#include "fvCFD.H"
#include "surfactantTransportModel.H"

namespace Foam
{

class CswDiff
:
    public surfactantTransportModel
{
    scalar rho_s_;
    scalar D_;
    dimensionedScalar Deff_;

    scalar rho_sw_;

public:
    TypeName("CswDiff");

    CswDiff(const dictionary& dict, foamAuxFields* aux);

    // CswDiff functions

    virtual void correct
    (
        const volScalarField& Sb,
        const surfaceScalarField& phib,
        const volScalarField& eps,
        const volScalarField& qb,
        const volScalarField& qs
    ) const;
};

} // End namespace Foam

#endif
