// READING TRANSPORTPROPERTIES DICT AND FLUID PROPERTIES

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

dimensionedScalar rho_a
(
    transportProperties.lookup("rho_a")
);

dimensionedScalar rho_b
(
    transportProperties.lookup("rho_b")
);

dimensionedScalar prst
(
    transportProperties.lookup("prst")
);

// READING INPUT FIELDS

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

volScalarField Sa
(
    IOobject
    (
        "Sa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);
Sa = scalar(1.0) - Sb;

Info<< "Reading field K\n" << endl;

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

Info<< "Reading field eps (porosity)\n" << endl;

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

Info<< "Reading field Ua\n" << endl;

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

Info<< "Reading field Ub\n" << endl;

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

// RELATIVE PERMEABILITY MODEL AND FIELDS

autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New(transportProperties);

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

// EVALUATE FLUXES AT CELL FACES

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

// FOAM MODEL FIELDS
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

    // Nca
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
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Fshear
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

    // Fsurf
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
    // nD field 
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

    // [rg - rc] field -- net foam texture generation (generation - coalescence)
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

    // reactive term
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
    *foamAux.F = -fvc::div(phia); // + qa // TODO

    // acumulation coefficient field
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

    // local-equilibrium foam texture field
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
}
else if (foamModelName == "MehrabiFoam")
{   
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
}

// SURFACTANT CONCENTRATION TRANSPORT MODEL

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
    // Fcsw = -fvc::div(phib); // + qb // TODO

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

// mobility and fractional flow update values
krModel->correct(kra, krb, Sb, *foamAux.Cs); 
surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));
surfaceScalarField krbf ("krbf",fvc::interpolate(krb,"krb"));

surfaceScalarField Maf ("Maf",Kf*kraf/mu_a);
surfaceScalarField Mbf ("Mbf",Kf*krbf/mu_b);
surfaceScalarField Mf ("Mf",Maf+Mbf);
surfaceScalarField Faf ("Faf",Maf/Mf);
surfaceScalarField Fbf ("Fbf",Mbf/Mf);

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
Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) );

Ub.correctBoundaryConditions();  
Ua.correctBoundaryConditions();

// total darcy velocity and total cell face flux
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
    Ua + Ub
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

// GRAVITATIONAL EFFECTS

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

// gravitational effect mobility part update values
surfaceScalarField Laf ("Laf",rho_a*Kf*kraf/mu_a);
surfaceScalarField Lbf ("Lbf",rho_b*Kf*krbf/mu_b);
surfaceScalarField Lf ("Lf",Laf+Lbf);

surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());

// CAPILLARY PRESSURE MODELS

Info<< "\nInitializing Pc" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);

Info<< "\nInitializing dpcds" << endl;
volScalarField dpcds
(
    IOobject
    (
        "dpcds",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);

// // Class capillary pressure model
autoPtr<capillaryPressureModel> capPressModel = capillaryPressureModel::New(transportProperties);
capPressModel->correct(pc, dpcds, Sb);

surfaceScalarField dpcdsf ("dpcdsf",fvc::interpolate(dpcds,"dpcds"));
surfaceScalarField phiPc("phiPc",Mbf * dpcdsf * fvc::snGrad(Sb) * mesh.magSf());

// SOURCE TERMS (WELL MODELS)

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

autoPtr<wellModel> wellModel = wellModel::New(wellProperties);
