Info<< "Reading field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Sa\n" << endl;

volScalarField Sa
(
    IOobject
    (
        "Sa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Sb\n" << endl;

volScalarField Sb
(
    IOobject
    (
        "Sb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Sc
(
    IOobject
    (
        "Sc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);
Sc = scalar(1.0) - Sa - Sb;

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar mu_a
(
    transportProperties.lookup("mu_a")
);

dimensionedScalar mu_b
(
    transportProperties.lookup("mu_b")
);

dimensionedScalar mu_c
(
    transportProperties.lookup("mu_c")
);

dimensionedScalar rho_a
(
    transportProperties.lookup("rho_a")
);

dimensionedScalar rho_b
(
    transportProperties.lookup("rho_b")
);

dimensionedScalar rho_c
(
    transportProperties.lookup("rho_c")
);

dimensionedScalar prst
(
    transportProperties.lookup("prst")
);

// Class relative permeability model

autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New(transportProperties);

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField Kf(fvc::interpolate(K,"K"));

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    prst
);

volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField Ub
(
    IOobject
    (
        "Ub",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField Uc
(
    IOobject
    (
        "Uc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phib
(
    IOobject
    (
        "phib",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(Ub)
);

surfaceScalarField phia
(
    IOobject
    (
        "phia",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(Ua)
);

surfaceScalarField phic
(
    IOobject
    (
        "phic",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(Uc)
);

volScalarField kra
(
    IOobject
    (
        "kra",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);

volScalarField krb
(
    IOobject
    (
        "krb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);

volScalarField krc
(
    IOobject
    (
        "krc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);

// CREATE zero RHS of equations 
dimensionedScalar zeroRHS
(
    "zero",
    Sb.dimensions()/dimTime,
    0.0
);

// FOAM RELATED FIELDS
foamAuxFields foamAux;

word foamModelName(transportProperties.lookup("foamModel")); 

if (foamModelName == "STARS")
{
    // Fdry
    foamAux.Fdry = new volScalarField
    (
        IOobject
        (
            "Fdry",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    // Fshear
    foamAux.Nca = new volScalarField
    (
        IOobject
        (
            "Nca",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 1.0)
    );

    foamAux.Fshear = new volScalarField 
    (
        IOobject
        (
            "Fshear",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    foamAux.Fsurf = new volScalarField 
    (
        IOobject
        (
            "Fsurf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    foamAux.Foil = new volScalarField 
    (
        IOobject
        (
            "Foil",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    // MRF
    foamAux.MRF = new volScalarField
    (
        IOobject
        (
            "MRF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

} else if (foamModelName == "ashoori")
{
    Info<< "Reading field nD\n" << endl;

    foamAux.nD = new volScalarField
    (
        IOobject
        (
            "nD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    foamAux.rgrc = new volScalarField
    (
        IOobject
        (
            "rgrc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        0.0
    );

    foamAux.F = new volScalarField
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    *foamAux.F = -fvc::div(phia); // + qa

    foamAux.AcumCoeff_nD = new volScalarField
    (
        IOobject
        (
            "AcumCoeff_nD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    *foamAux.AcumCoeff_nD = eps*(scalar(1.0) - Sb + VSMALL);

    foamAux.nDLE = new volScalarField
    (
        IOobject
        (
            "nDLE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    foamAux.MRF = new volScalarField
    (
        IOobject
        (
            "MRF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );
}
else if (foamModelName == "MehrabiFoam")
{
    foamAux.MRF = new volScalarField
    (
        IOobject
        (
            "MRF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );
}

// SURFACTANT CONCENTRATION TRANSPORT 

word surfTranspModelName(transportProperties.lookup("surfactantTransportModel")); 

if (surfTranspModelName == "Csw")
{
    Info<< "Reading field Cs\n" << endl;

    foamAux.Cs = new volScalarField
    (
        IOobject
        (
            "Cs",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    foamAux.Fcsw = new volScalarField
    (
        IOobject
        (
            "Fcsw",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    // Fcsw = -fvc::div(phib); // + qb

    foamAux.AcumCoeff_csw = new volScalarField
    (
        IOobject
        (
            "AcumCoeff_csw",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    // AcumCoeff_csw = eps*Sb + (scalar(1.0) - eps)*rho_sw*Kd; // Kd = dpcsds


    // Adsorption fields
    foamAux.CsEq = new volScalarField
    (
        IOobject
        (
            "CsEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    foamAux.dCsEqdCs = new volScalarField
    (
        IOobject
        (
            "dCsEqdCs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    foamAux.oneField = new volScalarField
    (
        IOobject
        (
            "oneField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );
} 

autoPtr<surfactantTransportModel> surfTranspModel = surfactantTransportModel::New(transportProperties, &foamAux);
autoPtr<foamModel> foamModel = foamModel::New(transportProperties, &foamAux);


// mobility and fractional flow 
krModel->correct(kra, krb, krc, Sa, Sb, *foamAux.Cs); 
surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));
surfaceScalarField krbf ("krbf",fvc::interpolate(krb,"krb"));
surfaceScalarField krcf ("krcf",fvc::interpolate(krc,"krc"));

surfaceScalarField Maf ("Maf",Kf*kraf/mu_a);
surfaceScalarField Mbf ("Mbf",Kf*krbf/mu_b);
surfaceScalarField Mcf ("Mcf",Kf*krcf/mu_c);
surfaceScalarField Mf ("Mf",Maf+Mbf+Mcf);
surfaceScalarField Faf ("Faf",Maf/Mf);
surfaceScalarField Fbf ("Fbf",Mbf/Mf);

volScalarField Fa
(
    IOobject
    (
        "Fa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);
Fa = (kra/mu_a) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );

volScalarField Fb
(
    IOobject
    (
        "Fb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);
Fa = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );

Ub.correctBoundaryConditions();  
Ua.correctBoundaryConditions();
Uc.correctBoundaryConditions();

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Ua + Ub + Uc
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

surfaceScalarField phiP = phi;

// ===========================

// GRAVITATIONAL EFFECTS
Info<< "\nReading g" << endl;

dimensionedVector g_vector
(
    transportProperties.lookup("g_vector")
);

surfaceVectorField g
(
    IOobject
    (
        "g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    g_vector
);

surfaceScalarField Laf ("Laf",rho_a*Kf*kraf/mu_a);
surfaceScalarField Lbf ("Lbf",rho_b*Kf*krbf/mu_b);
surfaceScalarField Lcf ("Lcf",rho_c*Kf*krcf/mu_c);
surfaceScalarField Lf ("Lf",Laf+Lbf+Lcf);

surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());

// CAPILLARY PRESSURE EFFECTS

// Info<< "\nInitializing Pc" << endl;
// volScalarField pc
// (
//     IOobject
//     (
//         "pc",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     0.0
// );

// Info<< "\nInitializing dpcds" << endl;
// volScalarField dpcds
// (
//     IOobject
//     (
//         "dpcds",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     0.0
// );

// // // Class capillary pressure model
// autoPtr<capillaryPressureModel> capPressModel = capillaryPressureModel::New(transportProperties);
// capPressModel->correct(pc, dpcds, Sb);

// surfaceScalarField dpcdsf ("dpcdsf",fvc::interpolate(dpcds,"dpcds"));
// surfaceScalarField phiPc("phiPc",Mbf * dpcdsf * fvc::snGrad(Sb) * mesh.magSf());


// SOURCE TERMS

// volume rate source term [m3 / s m3] -> [1 / s]
volScalarField qt
(
    IOobject
    (
        "qt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,0)
);

volScalarField qa
(
    IOobject
    (
        "qa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField qb
(
    IOobject
    (
        "qb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField qt_inj
(
    IOobject
    (
        "qt_inj",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Foam::mag(max(qt, scalar(0)))
);

volScalarField qt_prod
(
    IOobject
    (
        "qt_prod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Foam::mag(min(qt, scalar(0)))
);

Info<< "Reading wellProperties\n" << endl;

IOdictionary wellProperties
(
    IOobject
    (
        "wellProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Switch wellActive(wellProperties.headerOk());
// Info << "wellActive: " << wellActive << nl << endl;

autoPtr<wellModel> wellModel = wellModel::New(wellProperties);


foamAux.printStatus();