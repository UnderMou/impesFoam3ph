/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar Gd = 0.0;
scalar Epsilon = 0.0;
scalar meanGd = 0.0;
scalar meanEpsilon = 0.0;

if (mesh.nInternalFaces())
{
    scalarField Gd_values
        (   
            K*(rho_c.value()-rho_a.value())*mag(g) / (mag(U)*mu_c+VSMALL)
        );

    Gd = gMax(Gd_values);

    meanGd = 0.5*(gSum(Gd_values)/gSum(mesh.V().field()));

    // scalarField Pe_values
    //     (  
    //         mag(U)*mu_c.value()* 
    //         K*(rho_c.value()-rho_a.value())*mag(g) / (mag(U)*mu_c+VSMALL)
    //     );

    // scalarField Epsilon_values
    //     (   
    //         K*(rho_c.value()-rho_a.value())*mag(g) / (mag(U)*mu_c+VSMALL)
    //     );

    // Gd = gMax(Gd_values);

    // meanGd = 0.5*(gSum(Gd_values)/gSum(mesh.V().field()));
}

Info<< "Gd mean: " << meanGd << "," << " max: " << Gd << endl;
// Info<< "Epsilon mean: " << meanEpsilon << "," << " max: " << Epsilon << endl;

// CFLUse = CoNum;
// maxDeltaTFact = maxCo/(CFLUse + SMALL);

// ************************************************************************* //
