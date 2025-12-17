/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    impesFoam2ph

Description
    Solves two-phase flow in porous media accounting with foam and neglecting gravity
    and capillary effects

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "simpleControl.H"
#include "fvConstraints.H"

#include "RelPerm.H"
#include "FoamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    while (simple.loop(runTime))
    {

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        Info<< "\nCalculating p and Sb field\n" << endl;     

        while (simple.correctNonOrthogonal())
        {   

            correct_kra(kra,Sb,Sa_minValue,Sb_minValue,kra_maxValue,a_expValue);
            correct_krb(krb,Sb,Sa_minValue,Sb_minValue,krb_maxValue,b_expValue);

            // Foam model
            correct_Nca(Nca, U, mu_b, sigma_baValue); 
            correct_Fdry(Fdry, Sb, fmmobValue, SFValue, sfbetValue);
            correct_Fshear(Fshear, Nca, fmcapValue, epcapValue);
            correct_MRF(MRF, fmmobValue, Fdry, Fshear);
            correct_kraMRF(kra, MRF);

            kraf = fvc::interpolate(kra,"kra");
            krbf = fvc::interpolate(krb,"krb");
            Maf = Kf*kraf/mu_a;	
            Mbf = Kf*krbf/mu_b;
            Mf = Maf+Mbf;
            Fbf = Mbf/Mf;
            Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) );

            fvScalarMatrix pEqn
            (
                fvm::laplacian(-Mf, p)
                ==
                zeroRHS
            );

            pEqn.solve();

            phiP = pEqn.flux();

            phi = phiP;

            phib = Fbf*phiP;

            phia = phi - phib;

            U = fvc::reconstruct(phi);
            U.correctBoundaryConditions();

            Ub = fvc::reconstruct(phib);
            Ua = U-Ub;

            Ub.correctBoundaryConditions();  
            Ua.correctBoundaryConditions();

            forAll(mesh.boundary(),patchi)
            {
                if (isA< fixedValueFvPatchField<vector> >(Ua.boundaryField()[patchi]))
                {
                    phia.boundaryFieldRef()[patchi] = Ua.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
                }
                if (isA< fixedValueFvPatchField<vector> >(Ub.boundaryField()[patchi]))
                {
                    phib.boundaryFieldRef()[patchi] = Ub.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
                }
            }

            fvScalarMatrix SbEqn
            (
                eps*fvm::ddt(Sb) + fvc::div(phib) 
                ==
                zeroRHS
            );

            SbEqn.solve();

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
