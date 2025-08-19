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
#include "CapillaryPressure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    labelList prod_wells;
    forAll(qt, celli)
    {
        if (qt[celli] < 0)
        {
            prod_wells.append(celli);
        }
    }
    
    labelList inj_wells;
    forAll(qt, celli)
    {
        if (qt[celli] > 0)
        {
            inj_wells.append(celli);
        }
    }

    scalar zero = 1e-8;

    scalar minMf = gMin(Mf);
    surfaceScalarField Mf_safe = Mf + dimensionedScalar("SMALL", Mf.dimensions(), SMALL);
    scalar minLf = gMin(Lf);
    surfaceScalarField Lf_safe = Lf + dimensionedScalar("SMALL", Lf.dimensions(), SMALL);

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
        #include "GdEpsilon.H"

        // Info<< "\nCalculating p and Sb field\n" << endl;

        fvModels.correct();      

        while (simple.correctNonOrthogonal())
        {      
            // Pressure terms
            // Nca = correct_Nca(p, K, sigma_baValue);
            // Fshear = correct_Fshear(p, K, fmcapValue, epcapValue, sigma_baValue);
            // Foil = correct_Foil(Sa, Sb, Sa_minValue, Sb_minValue, fmoilValue, floilValue, epoilValue);
            kra = correct_kra(Sa,Sb,Sc,p,K,Sa_minValue,Sb_minValue,Sc_minValue,kra_maxValue,a_expValue,fmmobValue,SFValue,sfbetValue,U,fmcapValue,epcapValue,mu_b,sigma_baValue,fmoilValue,floilValue,epoilValue);
            krb = correct_krb(Sb,Sa_minValue,Sb_minValue,Sc_minValue,krb_maxValue,b_expValue);
            krc = correct_krc(Sc,Sa_minValue,Sb_minValue,Sc_minValue,krc_maxValue,c_expValue);
            kraf = fvc::interpolate(kra,"kra");
            krbf = fvc::interpolate(krb,"krb");
            krcf = fvc::interpolate(krc,"krc");
            Info<< "kri - OK\n" << endl;
            Maf = Kf*kraf/mu_a;	
            Mbf = Kf*krbf/mu_b;
            Mcf = Kf*krcf/mu_c;
            Mf = Maf+Mbf+Mcf;
            

            // Gravity terms
            Laf = rho_a*Kf*kraf/mu_a;
            Lbf = rho_b*Kf*krbf/mu_b;
            Lcf = rho_c*Kf*krcf/mu_c;
            Lf = Laf+Lbf+Lcf;
            // Info<< "Lif - OK\n" << endl;

            phiG = (Lf * g) & mesh.Sf();
            // Info<< "phiG - OK\n" << endl;

            minMf = gMin(Mf);
            Mf_safe = Mf + dimensionedScalar("SMALL", Mf.dimensions(), SMALL);
            Faf = Maf / Mf_safe;
            Fbf = Mbf / Mf_safe;

            // Info<< "Fif - OK\n" << endl;
            Fa = (kra/mu_a) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );
            Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );
            // Info<< "Fi - OK\n" << endl;
            Fc = 1.0 - Fa - Fb;
            Fcf = 1.0 - Faf - Fbf;

            // Capillary pressure terms
            pcgo = correct_pcgo(Sa,epsilon);
            pcwo = correct_pcwo(Sb,epsilon);

            Dpcwo_DSw = correct_Dpcwo_DSw(Sb,epsilon);
            Dpcgo_DSg = correct_Dpcgo_DSg(Sa,epsilon);
            Dpcwo_DSw_f = fvc::interpolate(Dpcwo_DSw,"Dpcwo_DSw");
            Dpcgo_DSg_f = fvc::interpolate(Dpcgo_DSg,"Dpcgo_DSg");

            phiPc_a = Maf * Dpcgo_DSg_f * fvc::snGrad(Sa) * mesh.magSf();
            phiPc_b = Mbf * Dpcwo_DSw_f * fvc::snGrad(Sb) * mesh.magSf();
            phiPc = phiPc_a + phiPc_b;
            // Info<< "phiPC - OK\n" << endl;
            
            forAll(mesh.boundary(),patchi)
            {   
                if ( Ua.boundaryField()[patchi].type() == "slip" )
                {   
                    phiG.boundaryFieldRef()[patchi] = 0.0;
                }
            }

            // Pressure equation
            fvScalarMatrix pEqn
            (
                fvm::laplacian(-Mf, p) + fvc::div(phiG) + fvc::div(phiPc)
                ==
                // fvModels.source(p)
                qt
            );
            // Info<< "pEqn - OK\n" << endl;

            pEqn.solve();

            phiP = pEqn.flux();

            forAll(mesh.boundary(),patchi)
            {   
                if ( Ua.boundaryField()[patchi].type() == "slip" )
                {   
                    phiP.boundaryFieldRef()[patchi] = 0.0;
                }
            }

            // Info<< "phiP - OK\n" << endl;

            phi = phiP + phiG + phiPc;

            minLf = gMin(Lf);
            Lf_safe = Lf + dimensionedScalar("SMALL", Lf.dimensions(), SMALL);

            phia = Faf*phiP + (Laf/Lf_safe)*phiG + phiPc_a;
            phib = Fbf*phiP + (Lbf/Lf_safe)*phiG + phiPc_b;
            phic = phi - phia - phib;

            // Info<< "phi_phases - OK\n" << endl;

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
                eps*fvm::ddt(Sa) + fvc::div(phia)
                ==
                fg_inj*qt_inj - Fa*qt_prod
                // fvModels.source(Sa)
            );

            SaEqn.solve(); 

            fvScalarMatrix SbEqn
            (
                eps*fvm::ddt(Sb) + fvc::div(phib)
                ==
                fw_inj*qt_inj - Fb*qt_prod
                // fvModels.source(Sb)
            );

            SbEqn.solve(); 

            Sa = correct_Sa(Sa, Sa_minValue, Sb_minValue, Sc_minValue);
            Info << "Saturation a: " << " Min(Sa) = " << gMin(Sa) << " Max(Sa) = " << gMax(Sa) << endl;

            Sb = correct_Sb(Sb, Sa_minValue, Sb_minValue, Sc_minValue);
            Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

            

            Sa.correctBoundaryConditions();  
            Sb.correctBoundaryConditions(); 

            Sc = oneField - Sa - Sb; 

            Sc = correct_Sc(Sc, Sa_minValue, Sb_minValue, Sc_minValue); 

            Info << "Saturation c: " << " Min(Sc) = " << gMin(Sc) << " Max(Sc) = " << gMax(Sc) << endl;

        }

        Info << runTime.writeTime() << endl;
        if (runTime.writeTime())
        {   
            // Info << "Foi!!!" << endl;
            std::ofstream outFile("production.txt", std::ios::app);
            forAll(prod_wells,i)
            {  
                outFile << runTime.timeName() << ","                                            // TimeStep
                << (Fb[inj_wells[i]]*qt_inj[inj_wells[i]]*mesh.V()[inj_wells[i]]) << ","        // qw_inj
                << (Fa[inj_wells[i]]*qt_inj[inj_wells[i]]*mesh.V()[inj_wells[i]]) << ","        // qg_inj
                << (Fb[prod_wells[i]]*qt_prod[prod_wells[i]]*mesh.V()[prod_wells[i]]) << ","    // qw_prod
                << (Fa[prod_wells[i]]*qt_prod[prod_wells[i]]*mesh.V()[prod_wells[i]]) << ","    // qg_prod
                << (Fc[prod_wells[i]]*qt_prod[prod_wells[i]]*mesh.V()[prod_wells[i]]) << "\n";  // qo_prod
            }
            outFile.close();
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
