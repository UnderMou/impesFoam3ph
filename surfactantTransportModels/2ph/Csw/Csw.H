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
    Foam::surfactantTransportModel::Csw

Description
    Surfactant transport by aqueous phase (only) model for two-phase flow.

    Reference:
    \verbatim
        Fritis, G. C., Paz, P. S., Lozano, L. F., & Chapiro, G. (2024). 
        On the Riemann problem for the foam displacement in porous media with linear adsorption. 
        SIAM Journal on Applied Mathematics, 84(2), 581-601.
    \endverbatim

    The default model coefficients are
    \verbatim

        Csw
        {

            rho_s 2.0e3;

            isothermModel Henry;

            Henry
            {
                Kd    0.05;
            }

        }

    \endverbatim

SourceFiles
    Csw.C

\*---------------------------------------------------------------------------*/

#ifndef Csw_H
#define Csw_H

#include "fvCFD.H"
#include "surfactantTransportModel.H"

namespace Foam
{

class Csw
:
    public surfactantTransportModel
{
    scalar rho_s_;

    scalar rho_sw_;

public:
    TypeName("Csw");

    Csw(const dictionary& dict, foamAuxFields* aux);

    // Csw functions

    virtual void correct
    (
        const volScalarField& Sb,
        const surfaceScalarField& phib,
        const volScalarField& eps
    ) const;
};

} // End namespace Foam

#endif
