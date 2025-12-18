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

// const dictionary& STARSDict = transportProperties.subDict("STARS");
// scalar fmmobValue = readScalar(STARSDict.lookup("fmmob"));
// scalar SFValue = readScalar(STARSDict.lookup("SF"));
// scalar sfbetValue = readScalar(STARSDict.lookup("sfbet"));
// scalar fmcapValue = readScalar(STARSDict.lookup("fmcap"));
// scalar epcapValue = readScalar(STARSDict.lookup("epcap"));
// scalar sigma_baValue = readScalar(STARSDict.lookup("sigma_ba"));

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

#include "createPhi.H"
surfaceScalarField phiP = phi;

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

krModel->correct(kra, krb, Sb);
surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));
surfaceScalarField krbf ("krbf",fvc::interpolate(krb,"krb"));

// mobility and fractional flow 
surfaceScalarField Maf ("Maf",Kf*kraf/mu_a);
surfaceScalarField Mbf ("Mbf",Kf*krbf/mu_b);
surfaceScalarField Mf ("Mf",Maf+Mbf);
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

// CREATE zero RHS of equations 
dimensionedScalar zeroRHS
(
    "zero",
    Sb.dimensions()/dimTime,
    0.0
);

// Class foam model

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
        dimensionedScalar("zero", dimless, 0.0)
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


autoPtr<foamModel> foamModel = foamModel::New(transportProperties, &foamAux);

// ===========================

// foamModel->correct(kra, U, Sb);
