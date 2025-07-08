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
            //- mobility computation 
            Nca = correct_Nca(p, K, sigma_baValue);
            Fshear = correct_Fshear(p, K, fmcapValue, epcapValue, sigma_baValue);
            Foil = correct_Foil(Sa, Sb, Sa_minValue, Sb_minValue, fmoilValue, floilValue, epoilValue);
            kra = correct_kra(Sa,Sb,p,K,Sa_minValue,Sb_minValue,Sc_minValue,kra_maxValue,a_expValue,fmmobValue,SFValue,sfbetValue,U,fmcapValue,epcapValue,mu_b,sigma_baValue,fmoilValue,floilValue,epoilValue);
            krb = correct_krb(Sb,Sa_minValue,Sb_minValue,Sc_minValue,krb_maxValue,b_expValue);
            krc = correct_krc(Sa,Sb,Sa_minValue,Sb_minValue,Sc_minValue,krc_maxValue,c_expValue);
            kraf = fvc::interpolate(kra,"kra");
            krbf = fvc::interpolate(krb,"krb");
            krcf = fvc::interpolate(krc,"krc");
            Maf = Kf*kraf/mu_a;	
            Mbf = Kf*krbf/mu_b;
            Mcf = Kf*krcf/mu_c;
            Mf = Maf+Mbf+Mcf;
            Faf = Maf/Mf;
            Fbf = Mbf/Mf;
            Fa = (kra/mu_a) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );
            Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );

            fvScalarMatrix pEqn
            (
                fvm::laplacian(-Mf, p)
                ==
                fvModels.source(p)
            );

            pEqn.solve();

            phiP = pEqn.flux();

            phi = phiP;

            phia = Faf*phiP;
            phib = Fbf*phiP;

            phic = phi - phia - phib;

            U = fvc::reconstruct(phi);
            U.correctBoundaryConditions();

            Ua = fvc::reconstruct(phia);
            Ub = fvc::reconstruct(phib);
            Uc = U-Ua-Ub;

            Ua.correctBoundaryConditions();  
            Ub.correctBoundaryConditions();
            Uc.correctBoundaryConditions();

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
                if (isA< fixedValueFvPatchField<vector> >(Uc.boundaryField()[patchi]))
                {
                    phic.boundaryFieldRef()[patchi] = Uc.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
                }
            }

            fvScalarMatrix SaEqn
            (
                eps*fvm::ddt(Sa) + fvm::div(phia, Sa) // fvc::div(phia) 
                ==
                fvModels.source(Sa)
            );

            SaEqn.solve(); 

            fvScalarMatrix SbEqn
            (
                eps*fvm::ddt(Sb) + fvm::div(phib, Sb) // fvc::div(phib) 
                ==
                fvModels.source(Sb)
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
