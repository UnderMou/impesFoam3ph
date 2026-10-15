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
    Solves two-phase flow in porous media through darcy law (Gravity and capillary effects
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

        #include "CourantNo.H"  // CFL for each phase
        #include "GdEpsilon.H"  // Gravity over convective effects
    
        while (simple.correctNonOrthogonal())
        {   
            
            // Relative permeability model
            Info<< "Using relative permeability model: " << krModel->type() << nl << endl;
            krModel->correct(kra, krb, Sb, *foamAux.Cs);

            // Foam model
            Info<< "Using foam model: " << foamModel->type() << nl << endl;
            foamModel->correct(kra, U, Sa, Sb, phia, eps, K, fvc::grad(p));
            Info<<"Foam ok"<<nl<<endl;
            // Total mobility (M_i) and fractional flux
            kraf = fvc::interpolate(kra,"kra");
            krbf = fvc::interpolate(krb,"krb");
            Maf = Kf*kraf/mu_a;	
            Mbf = Kf*krbf/mu_b;
            Mf = Maf+Mbf;

            Faf = Maf/Mf;
            Fbf = Mbf/Mf;
            Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) );
            mob_t = kra/mu_a + krb/mu_b;
            mob_a = kra/mu_a;
            mob_b = krb/mu_b;

            // Gravitational effects (L_i)
            Laf = rho_a*Kf*kraf/mu_a;
            Lbf = rho_b*Kf*krbf/mu_b;	
            Lf = Laf+Lbf;
            phiG = (Lf * g) & mesh.Sf();

            // Capillary pressure model
            Info<< "Using capillary pressure model: " << capPressModel->type() << nl << endl;
            capPressModel->correct(pc, dpcds, Sb);

            dpcdsf = fvc::interpolate(dpcds,"dpcds");
            // The BC must be updated before phiPc is evaluated because, when capillarity
            // is active, it determines snGrad(Sb).
            Sb.correctBoundaryConditions();
            phiPc = Mbf * dpcdsf * fvc::snGrad(Sb) * mesh.magSf();

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
            wellModel->source_pEqn(pEqn,p,mob_t,WI,wellCoeff,wellSource,rho_a.value(),rho_b.value(),mob_a,mob_b,g_vector,qt,qb);
            if (usePressureReference)
            {
                pEqn.setReference(pRefCell, pRefValue);
            }
            pEqn.solve();
            phiP = pEqn.flux();

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
                                phiP.boundaryField()[patchi] + phiG.boundaryField()[patchi] + phiPc.boundaryField()[patchi]
                            )
                        )
                        << endl;
                }
            }

            // phase fluxe at cell faces
            phib = Fbf*phiP + (Lbf/Lf)*phiG + phiPc;

            if(capPressModel->type() == "noCapillaryPressure")
            {
                forAll(mesh.boundary(), patchi)
                {
                    if (isA<wallFvPatch>(mesh.boundary()[patchi]))
                    {
                        phib.boundaryFieldRef()[patchi] = 0.0;
                    }
                }
            }
            phia = phi - phib;
            
            // In capillary-active cases this should already be ~0 from the BC-derived
            // gradient. In no-capillary/degenerate cases, the saturation gradient cannot
            // influence phib, so the zero-flux condition must be imposed directly.
            // forAll(mesh.boundary(), patchi)
            // {
            //     if
            //     (
            //         Sb.boundaryField()[patchi].type()
            //     == "darcyNoFluxSaturation"
            //     )
            //     {
            //         phib.boundaryFieldRef()[patchi] = 0.0;
            //     }
            // }



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
                }
            }

            

            // // total flux equals zero on walls for each phase
            // forAll(mesh.boundary(),patchi)
            // {   
            //     if ( Ua.boundaryField()[patchi].type() == "slip" )
            //     {   
            //         phia.boundaryFieldRef()[patchi] = 0.0;
            //         phib.boundaryFieldRef()[patchi] = 0.0;
            //     }
            // }

            // correct darcy velocities at boundaries
            Ua = fvc::reconstruct(phia);
            Ub = fvc::reconstruct(phib);
            Ub.correctBoundaryConditions();  
            Ua.correctBoundaryConditions();
            U = Ua + Ub; // Correct U according with Ua and Ub Boundary values
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

            // well model correction
            Info<< "Using well model: " << wellModel->type() << nl << endl;
            wellModel->correct(qt,qb,Fb,p,runTime.timeOutputValue(),mob_t,WI,p_bh,qs,*foamAux.Cs,rho_a.value(),rho_b.value(),mob_a,mob_b,g_vector);
            Info << "BHP = " << gMax(p_bh.internalField()) << endl;
            
            if (capPressModel->type() == "noCapillaryPressure")
            {
                forAll(mesh.boundary(), patchi)
                {
                    if (isA<wallFvPatch>(mesh.boundary()[patchi]))
                    {
                        phib.boundaryFieldRef()[patchi] = 0.0;
                    }
                }
            }
            forAll(mesh.boundary(), patchi)
            {
                if (isA<wallFvPatch>(mesh.boundary()[patchi]))
                {
                    Info<< "BEFORE SbEqn - "
                        << mesh.boundary()[patchi].name()
                        << " max|phib| = "
                        << gMax
                        (
                            mag(phib.boundaryField()[patchi])
                        )
                        << endl;
                }
            }
            // phase saturation equation
            fvScalarMatrix SbEqn
            (
                eps*fvm::ddt(Sb) + fvc::div(phib)
            );
            wellModel->source_SbEqn(SbEqn,Sb,Fb,p,runTime.timeOutputValue(),qb);
            // SbEqn.solve();

            // Sa = scalar(1.0) - Sb;

            // Sb.correctBoundaryConditions();

            // SbEqn.solve();
            const label debugCell = 5; // coloque aqui uma owner problemática

            const cell& cFaces = mesh.cells()[debugCell];

            forAll(cFaces, i)
            {
                const label facei = cFaces[i];

                if (facei < mesh.nInternalFaces())
                {
                    const label owner = mesh.faceOwner()[facei];
                    const label nei   = mesh.faceNeighbour()[facei];

                    label otherCell;

                    scalar flux;

                    if (owner == debugCell)
                    {
                        otherCell = nei;
                        flux = phib[facei];
                    }
                    else
                    {
                        otherCell = owner;
                        flux = -phib[facei];
                    }

                    const vector Cf = mesh.Cf()[facei];
                    const vector Sf = mesh.Sf()[facei];
                    const vector n  = Sf/mag(Sf);

                    Info<< nl
                        << "face        = " << facei << nl
                        << "Cf          = " << Cf << nl
                        << "Sf          = " << Sf << nl
                        << "normal      = " << n << nl
                        << "flux out    = " << flux << nl
                        << "other cell  = " << otherCell << nl
                        << "C neighbour = " << mesh.C()[otherCell] << nl
                        << "Sb owner    = " << Sb[debugCell] << nl
                        << "Sb neighbour= " << Sb[otherCell] << nl;
                }
            }

            volScalarField divPhibDebug
            (
                IOobject
                (
                    "divPhibDebug",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::div(phib)
            );

            Info<< "===== BEFORE SbEqn =====" << nl
                << "cell       = " << debugCell << nl
                << "C          = " << mesh.C()[debugCell] << nl
                << "Sb oldTime = " << Sb.oldTime()[debugCell] << nl
                << "Sb current = " << Sb[debugCell] << nl
                << "div(phib)  = " << divPhibDebug[debugCell] << nl
                << "eps        = " << eps[debugCell] << nl
                << "deltaT     = " << runTime.deltaTValue()
                << endl;

            scalar SbPred =
                Sb.oldTime()[debugCell]
            - runTime.deltaTValue()
            *divPhibDebug[debugCell]
            /eps[debugCell];

            Info<< "Sb predicted = " << SbPred << endl;

            SbEqn.solve();

            Info<< "===== AFTER SbEqn =====" << nl
                << "Sb solved = " << Sb[debugCell]
                << endl;

            forAll(mesh.boundary(), patchi)
            {
                if
                (
                    Sb.boundaryField()[patchi].type()
                    == "darcyNoFluxSaturation"
                )
                {
                    const fvPatchScalarField& Sbp =
                        Sb.boundaryField()[patchi];

                    scalarField SbOwner(Sbp.patchInternalField());

                    Info<< "BEFORE correctBoundaryConditions - "
                        << mesh.boundary()[patchi].name()
                        << nl
                        << "max|Sb_face - Sb_owner| = "
                        << gMax(mag(Sbp - SbOwner))
                        << endl;
                }
            }

            // const cell& cFaces = mesh.cells()[debugCell];

            scalar sumPhi = 0.0;

            Info<< "DEBUG cell = " << debugCell
                << " nFaces = " << cFaces.size()
                << endl;

            forAll(cFaces, i)
            {
                const label facei = cFaces[i];

                scalar flux = 0.0;

                if (facei < mesh.nInternalFaces())
                {
                    // Internal face
                    if (mesh.faceOwner()[facei] == debugCell)
                    {
                        flux = phib[facei];
                    }
                    else
                    {
                        flux = -phib[facei];
                    }

                    Info<< "internal face = " << facei
                        << " flux(outward) = " << flux
                        << endl;
                }
                else
                {
                    // Boundary face
                    const label patchi =
                        mesh.boundaryMesh().whichPatch(facei);

                    Info<< "boundary face = " << facei
                        << " patchi = " << patchi
                        << endl;

                    if (patchi < 0 || patchi >= mesh.boundary().size())
                    {
                        Info<< "ERROR: invalid patch for face "
                            << facei << endl;

                        continue;
                    }

                    const fvPatch& patch =
                        mesh.boundary()[patchi];

                    const label localFace =
                        facei - patch.start();

                    Info<< "    patch = " << patch.name()
                        << " localFace = " << localFace
                        << " patchSize = " << patch.size()
                        << endl;

                    if (localFace < 0 || localFace >= patch.size())
                    {
                        Info<< "ERROR: invalid localFace" << endl;
                        continue;
                    }

                    flux =
                        phib.boundaryField()[patchi][localFace];

                    Info<< "    flux(outward) = "
                        << flux << endl;
                }

                sumPhi += flux;
            }

            Info<< "sumPhi = " << sumPhi << endl;


            Sb.correctBoundaryConditions();

            // const label debugCell = 5;

            Info<< nl
                << "==============================" << nl
                << "BOUNDARIES OF CELL " << debugCell << nl
                << "Cell centre = " << mesh.C()[debugCell] << nl
                << "Sb owner    = " << Sb[debugCell] << nl;

            forAll(mesh.boundary(), patchi)
            {
                const fvPatch& patch = mesh.boundary()[patchi];

                const labelUList& faceCells = patch.faceCells();

                forAll(faceCells, facei)
                {
                    if (faceCells[facei] == debugCell)
                    {
                        Info<< "patch      = " << patch.name() << nl
                            << "face        = " << facei << nl
                            << "Cf          = "
                            << mesh.Cf().boundaryField()[patchi][facei] << nl
                            << "Sb face     = "
                            << Sb.boundaryField()[patchi][facei] << nl
                            << "Sb owner    = "
                            << Sb[debugCell] << nl
                            << "difference  = "
                            << Sb.boundaryField()[patchi][facei]
                            - Sb[debugCell]
                            << nl;
                    }
                }

            }

            Info<< "==============================" << endl;


            forAll(mesh.boundary(), patchi)
            {
                if
                (
                    Sb.boundaryField()[patchi].type()
                    == "darcyNoFluxSaturation"
                )
                {
                    const fvPatchScalarField& Sbp =
                        Sb.boundaryField()[patchi];

                    scalarField SbOwner(Sbp.patchInternalField());

                    Info<< "AFTER correctBoundaryConditions - "
                        << mesh.boundary()[patchi].name()
                        << nl
                        << "max|Sb_face - Sb_owner| = "
                        << gMax(mag(Sbp - SbOwner))
                        << endl;
                }
            }

            Sa = scalar(1.0) - Sb;
            Sa.correctBoundaryConditions();

            forAll(mesh.boundary(), patchi)
            {
                if
                (
                    Sb.boundaryField()[patchi].type()
                    == "darcyNoFluxSaturation"
                )
                {
                    const fvPatchScalarField& Sbp =
                        Sb.boundaryField()[patchi];

                    scalarField gradSb(Sbp.snGrad());

                    Info<< mesh.boundary()[patchi].name()
                        << " max|snGrad(Sb)| = "
                        << gMax(mag(gradSb))
                        << endl;
                }
            }

            Info << "Saturation a: " << " Min(Sa) = " << gMin(Sa) << " Max(Sa) = " << gMax(Sa) << endl;
            Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

            // Surfactant transport model
            Info<< "Using surfactant transport model: " << surfTranspModel->type() << nl << endl;
            surfTranspModel->correct(Sb, phib, eps, qb, qs);

        }

        Sb.correctBoundaryConditions();

        Sa = scalar(1.0) - Sb;
        Sa.correctBoundaryConditions();

        

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
