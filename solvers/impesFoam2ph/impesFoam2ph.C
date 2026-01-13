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
#include "wallFvPatch.H"

#include "relativePermeabilityModel.H"
#include "capillaryPressureModel.H"
#include "foamModel.H"
#include "surfactantTransportModel.H"
#include "isothermModel.H"
#include "wellModel.H"

// #include "wellModeling.H"

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
            // Relative permeability model
            Info<< "Using relative permeability model: " << krModel->type() << nl << endl;
            krModel->correct(kra, krb, Sb);

            // Foam model
            Info<< "Using foam model: " << foamModel->type() << nl << endl;
            foamModel->correct(kra, U, Sa, Sb, phia, eps);

            kraf = fvc::interpolate(kra,"kra");
            krbf = fvc::interpolate(krb,"krb");
            Maf = Kf*kraf/mu_a;	
            Mbf = Kf*krbf/mu_b;
            Mf = Maf+Mbf;

            Laf = rho_a*Kf*kraf/mu_a;
            Lbf = rho_b*Kf*krbf/mu_b;	
            Lf = Laf+Lbf;
            phiG = (Lf * g) & mesh.Sf();

            Fbf = Mbf/Mf;
            Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) );
            
            // Capillary pressure model
            Info<< "Using capillary pressure model: " << capPressModel->type() << nl << endl;
            capPressModel->correct(pc, dpcds, Sb);

            dpcdsf = fvc::interpolate(dpcds,"dpcds");
            phiPc = Mbf * dpcdsf * fvc::snGrad(Sb) * mesh.magSf();

            forAll(mesh.boundary(),patchi)
            {   
                if ( Ua.boundaryField()[patchi].type() == "slip" )
                {   
                    phiG.boundaryFieldRef()[patchi] = 0.0;
                }
            }

            fvScalarMatrix pEqn
            (
                fvm::laplacian(-Mf, p) + fvc::div(phiG) + fvc::div(phiPc)
                ==
                // zeroRHS
                qt
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

            phi = phiP + phiG + phiPc;

            phib = Fbf*phiP + (Lbf/Lf)*phiG + phiPc;

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

            // well correction
            wellModel->correct(qb,Fb,qt_inj,qt_prod);

            fvScalarMatrix SbEqn
            (
                eps*fvm::ddt(Sb) + fvc::div(phib) 
                ==
                // zeroRHS
                qb
            );

            SbEqn.solve();

            Sa = scalar(1.0) - Sb;

            Info << "Saturation a: " << " Min(Sa) = " << gMin(Sa) << " Max(Sa) = " << gMax(Sa) << endl;
            Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

            // Surfactant transport model
            // Info<< "Using isotherm model: " << isothermModel->type() << nl << endl;
            // isothermModel->correct();
            Info<< "Using surfactant transport model: " << surfTranspModel->type() << nl << endl;
            surfTranspModel->correct(Sb, phib, eps);

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
