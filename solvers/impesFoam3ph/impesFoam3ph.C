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
    impesFoam3ph

Description
    Solves three-phase flow in porous media through darcy law (Gravity and capillary effects
    included) incorporating foam effects (implicit and explicit texture), surfactant transport 
    and adsorption, and well models. 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "simpleControl.H"
#include "fvConstraints.H"
#include "wallFvPatch.H"

#include "relativePermeabilityModel.H"
#include "capillaryPressureModel.H"
#include "foamModel.H"
#include "surfactantTransportModel.H"
#include "isothermModel.H"
#include "wellModel.H"

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

        // Well inputs at this time step
        wellModel->updateWellInputs(runTime.value());

        #include "CourantNo.H"  // CFL for each phase
        #include "GdEpsilon.H"  // Gravity over convective effects
    
        while (simple.correctNonOrthogonal())
        {   
            // Relative permeability model
            Info<< "Using relative permeability model: " << krModel->type() << nl << endl;
            krModel->correct(kra, krb, krc, Sa, Sb, *foamAux.Cs);

            // Foam model
            Info<< "Using foam model: " << foamModel->type() << nl << endl;
            foamModel->correct(kra, U, Sa, Sb, Sc, phia, eps, K, fvc::grad(p));

            // Total mobility  (M_i) and fractional flux
            kraf = fvc::interpolate(kra,"kra");
            krbf = fvc::interpolate(krb,"krb");
            krcf = fvc::interpolate(krc,"krc");
            Maf = Kf*kraf/mu_a;	
            Mbf = Kf*krbf/mu_b;
            Mcf = Kf*krcf/mu_c;
            Mf = Maf+Mbf+Mcf;
            Mf = max(Mf, VSMALL);

            Faf = Maf/Mf;
            Fbf = Mbf/Mf;
            Fa = (kra/mu_a) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );
            Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );
            mob_t = kra/mu_a + krb/mu_b + krc/mu_c;
            mob_a = kra/mu_a;
            mob_b = krb/mu_b;
            mob_c = krc/mu_c;

            // Gravitational effects (L_i)
            Laf = rho_a*Kf*kraf/mu_a;
            Lbf = rho_b*Kf*krbf/mu_b;
            Lcf = rho_c*Kf*krcf/mu_c;
            Lf = Laf+Lbf+Lcf;
            Lf = max(Lf, VSMALL);
            phiG = (Lf * g) & mesh.Sf();

            // Capillary pressure model
            Info<< "Using capillary pressure model: " << capPressModel->type() << nl << endl;
            capPressModel->correct(pcow, pcgo,
                                   dpcow_dsw, dpcow_dsg,
                                   dpcgo_dsw, dpcgo_dsg,
                                   Sa, Sb, zeroField);

            dpcow_dsw_f = fvc::interpolate(dpcow_dsw,"dpcow_dsw");
            dpcow_dsg_f = fvc::interpolate(dpcow_dsg,"dpcow_dsg");
            dpcgo_dsw_f = fvc::interpolate(dpcgo_dsw,"dpcgo_dsw");
            dpcgo_dsg_f = fvc::interpolate(dpcgo_dsg,"dpcgo_dsg");
            
            Sa.correctBoundaryConditions();
            Sb.correctBoundaryConditions();

            phiPc_a = Maf * (dpcgo_dsw_f * fvc::snGrad(Sb) * mesh.magSf() + dpcgo_dsg_f * fvc::snGrad(Sa) * mesh.magSf());
            phiPc_b = Mbf * (dpcow_dsw_f * fvc::snGrad(Sb) * mesh.magSf() + dpcow_dsg_f * fvc::snGrad(Sa) * mesh.magSf());
            phiPc = phiPc_a - phiPc_b;

            // // zero gravitational on walls
            // forAll(mesh.boundary(),patchi)
            // {   
            //     if ( Ua.boundaryField()[patchi].type() == "slip" )
            //     {   
            //         phiG.boundaryFieldRef()[patchi] = 0.0;
            //     }
            // }

            // pressure equation
            fvScalarMatrix pEqn
            (
                fvm::laplacian(-Mf, p) + fvc::div(phiG) + fvc::div(phiPc)
            );
            wellModel->source_pEqn(pEqn,p,mob_t,WI,wellCoeff,wellSource,rho_a.value(),rho_b.value(),rho_c.value(),mob_a,mob_b,mob_c,g_vector,qt,qa,qb);
            if (usePressureReference)
            {
                pEqn.setReference(pRefCell, pRefValue);
            }
            pEqn.solve();
            phiP = pEqn.flux();

            // forAll(mesh.boundary(),patchi)
            // {   
            //     if ( Ua.boundaryField()[patchi].type() == "slip" )
            //     {   
            //         phiP.boundaryFieldRef()[patchi] = 0.0;
            //     }
            // }

            // total flux at cell faces
            phi = phiP + phiG + phiPc;

            // CHECK phi = 0 walls
            forAll(mesh.boundary(), patchi)
            {
                if
                (
                    p.boundaryField()[patchi].type()
                    == "darcyNoFluxPressure"
                )
                {
                    Info<< mesh.boundary()[patchi].name()
                        << " max|phiP + phiG + phiPc| = "
                        << gMax
                        (
                            mag
                            (
                                phiP.boundaryField()[patchi]
                            + phiG.boundaryField()[patchi]
                            + phiPc.boundaryField()[patchi]
                            )
                        )
                        << endl;
                }
            }

            // phase fluxe at cell faces
            phia = Faf*phiP + (Laf/Lf)*phiG + phiPc_a;
            phib = Fbf*phiP + (Lbf/Lf)*phiG - phiPc_b;
            
            if(capPressModel->type() == "noCapillaryPressure")
            {
                forAll(mesh.boundary(), patchi)
                {
                    if (isA<wallFvPatch>(mesh.boundary()[patchi]))
                    {
                        phia.boundaryFieldRef()[patchi] = 0.0;
                        phib.boundaryFieldRef()[patchi] = 0.0;
                    }
                }
            }
            phic = phi - phia - phib;

            forAll(mesh.boundary(), patchi)
            {
                if
                (
                    p.boundaryField()[patchi].type()
                    == "darcyNoFluxPressure"
                )
                {
                    Info<< mesh.boundary()[patchi].name()
                        << " max|phib| = "
                        << gMax
                        (
                            mag
                            (
                                phib.boundaryField()[patchi]
                            )
                        )
                        << endl;
                    Info<< mesh.boundary()[patchi].name()
                        << " max|phia| = "
                        << gMax
                        (
                            mag
                            (
                                phia.boundaryField()[patchi]
                            )
                        )
                        << endl;
                    Info<< mesh.boundary()[patchi].name()
                        << " max|phic| = "
                        << gMax
                        (
                            mag
                            (
                                phic.boundaryField()[patchi]
                            )
                        )
                        << endl;
                }
            }

            // correct darcy velocities at boundaries
            // U = fvc::reconstruct(phi);
            // U.correctBoundaryConditions();
            Ua = fvc::reconstruct(phia);
            Ub = fvc::reconstruct(phib);
            Uc = fvc::reconstruct(phic);
            // Uc = U-Ua-Ub;
            Ub.correctBoundaryConditions();  
            Ua.correctBoundaryConditions();
            Uc.correctBoundaryConditions();
            U = Ua + Ub + Uc; // Correct U according with Ua , Ub and Uc Boundary values
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

            // well model correction
            Info<< "Using well model: " << wellModel->type() << nl << endl;
            wellModel->correct(qt,qa,qb,Fa,Fb,p,runTime.timeOutputValue(),mob_t,WI,p_bh,qs,*foamAux.Cs,rho_a.value(),rho_b.value(),rho_c.value(),mob_a,mob_b,mob_c,g_vector);

            // phase saturation equation
            fvScalarMatrix SaEqn
            (
                eps*fvm::ddt(Sa) + fvc::div(phia)
            );
            wellModel->source_SaEqn(SaEqn,Sa,Fa,p,runTime.timeOutputValue(),qa);
            SaEqn.solve();

            fvScalarMatrix SbEqn
            (
                eps*fvm::ddt(Sb) + fvc::div(phib) 
            );
            wellModel->source_SbEqn(SbEqn,Sb,Fb,p,runTime.timeOutputValue(),qb);
            SbEqn.solve();

            // Sa = max(min(Sa, krModel.Sa_max_), krModel.Sa_min_);
            // Sb = max(min(Sb, krModel.Sb_max_), krModel.Sb_min_);
            Sa.correctBoundaryConditions();
            Sb.correctBoundaryConditions();
            Sc = scalar(1.0) - Sa - Sb;

            // Sc = min(max(Sc, scalar(0)), scalar(1));


            Info << "Saturation a: " << " Min(Sa) = " << gMin(Sa) << " Max(Sa) = " << gMax(Sa) << endl;
            Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;
            Info << "Saturation c: " << " Min(Sc) = " << gMin(Sc) << " Max(Sc) = " << gMax(Sc) << endl;

            // Surfactant transport model
            Info<< "Using surfactant transport model: " << surfTranspModel->type() << nl << endl;
            surfTranspModel->correct(Sb, phib, eps, qb, qs);

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
