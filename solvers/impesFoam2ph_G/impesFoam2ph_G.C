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
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "simpleControl.H"
#include "fvConstraints.H"
// #include "uniformDimensionedVectorField.H"
// #include "dimensionedVector.H"
#include "wallFvPatch.H"

#include "RelPerm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

   

    #include "createFields.H"

   

    // dimensionedScalar adm_1 = dimensionedScalar("adm_1",dimLength/dimTime,1.0);
    // dimensionedScalar adm_2 = dimensionedScalar("adm_2",dimLength*dimLength/dimTime,1.0);
    // dimensionedScalar adm_3 = dimensionedScalar("adm_2",dimTime,1.0);

    // // INITIALIZE K
    // forAll(K, i)
    // {   
    //     const scalar x = mesh.C()[i][0];
    //     const scalar y = mesh.C()[i][1];

    //     // K[i] = 1.0e-12 + 1.0e-13*x - 1.0e-13*y;
    //     // Info << i << "," << K[i] << "," << x << "," << y << endl;

    //     Info << i << "," << x << "," << y << endl;
    // }  

    while (simple.loop(runTime))
    {

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        Info<< "\nCalculating p and Sb field\n" << endl;

        fvModels.correct();      

        while (simple.correctNonOrthogonal())
        {      
            // Sb = correct_Sb(Sb, Sa_minValue, Sb_minValue); 
            // forAll(Sb, i)
            // {
            //     Info << Sb[i] << " ";
            // }
            // Info << endl;
            // Sb = correct_Sb(Sb, Sa_minValue, Sb_minValue);
            // Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;
            
            kra = correct_kra(Sb,Sa_minValue,Sb_minValue,kra_maxValue,a_expValue,fmmobValue,SFValue,sfbetValue,U,fmcapValue,epcapValue,mu_b,sigma_baValue);
            krb = correct_krb(Sb,Sa_minValue,Sb_minValue,krb_maxValue,b_expValue);
            kraf = fvc::interpolate(kra,"kra");
            krbf = fvc::interpolate(krb,"krb");
            Maf = Kf*kraf/mu_a;	
            Laf = rho_a*Kf*kraf/mu_a;	
            Mbf = Kf*krbf/mu_b;
            Lbf = rho_b*Kf*krbf/mu_b;
            Mf = Maf+Mbf;
            Lf = Laf+Lbf;

            phiG = (Lf * g) & mesh.Sf();

            Fbf = Mbf/Mf;
            Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) );

            forAll(mesh.boundary(),patchi)
            {   
                if ( Ua.boundaryField()[patchi].type() == "slip" )
                {   
                    phiG.boundaryFieldRef()[patchi] = 0.0;
                    // Info << Ua.boundaryField()[patchi].name() << endl;
                }
            }

            fvScalarMatrix pEqn
            (
                fvm::laplacian(-Mf, p) + fvc::div(phiG)
                ==
                fvModels.source(p) 
            );

            pEqn.solve();

            phiP = pEqn.flux();

            forAll(mesh.boundary(),patchi)
            {   
                if ( Ua.boundaryField()[patchi].type() == "slip" )
                {   
                    phiP.boundaryFieldRef()[patchi] = 0.0;
                }
            }

            phi = phiP + phiG;

            phib = Fbf*phiP + (Lbf/Lf)*phiG;

            phia = phi - phib;

            // int cont = 0;
            // forAll(mesh.boundary(),patchi)
            // {   
            //     //Info << Ua.boundaryField()[patchi] << endl;
            //     if ( Ua.boundaryField()[patchi].type() == "slip" )
            //     {   
            //         cont++;
            //         phia.boundaryFieldRef()[patchi] = 0.0;
            //         phi.boundaryFieldRef()[patchi] = 0.0;
            //     }
            //     if ( Ub.boundaryField()[patchi].type() == "slip" )
            //     {   
            //         cont++;
            //         phib.boundaryFieldRef()[patchi] = 0.0;
            //     }

            //     // cont++;
            // }
            // Info << cont << endl;

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
                fvModels.source(Sb)
            );

            SbEqn.solve();

            Sb = correct_Sb(Sb, Sa_minValue, Sb_minValue);
            Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;
            
            
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
