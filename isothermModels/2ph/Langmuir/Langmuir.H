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
    Foam::isothermModel::Langmuir

Description
    Langmuir adsorption model for two-phase flow.

    Reference:
    \verbatim
        Kalam, S., Abu-Khamsin, S. A., Kamal, M. S., & Patil, S. (2021).
        Surfactant adsorption isotherms: A review.
        ACS omega, 6(48), 32342-32348.
    \endverbatim

    The default model coefficients are
    \verbatim

        Langmuir
        {
            qo    0.3;
            KL    1.0;
        }

    \endverbatim

SourceFiles
    Langmuir.C

\*---------------------------------------------------------------------------*/

#ifndef Langmuir_H
#define Langmuir_H

#include "isothermModel.H"

namespace Foam
{

class Langmuir
:
    public isothermModel
{
    scalar qo_;
    scalar KL_;

public:
    TypeName("Langmuir");

    Langmuir(const dictionary& dict, foamAuxFields* aux);

    virtual void correct
    (
       
    ) const;
};

} // End namespace Foam

#endif
