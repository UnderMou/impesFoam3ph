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

dimensionedScalar Sa_min = transportProperties.lookup("Sa_min");
scalar Sa_minValue = Sa_min.value();

dimensionedScalar Sb_min = transportProperties.lookup("Sb_min");
scalar Sb_minValue = Sb_min.value();

dimensionedScalar a_exp = transportProperties.lookup("a_exp");
scalar a_expValue = a_exp.value();

dimensionedScalar b_exp = transportProperties.lookup("b_exp");
scalar b_expValue = b_exp.value();

dimensionedScalar kra_max = transportProperties.lookup("kra_max");
scalar kra_maxValue = kra_max.value();

dimensionedScalar krb_max = transportProperties.lookup("krb_max");
scalar krb_maxValue = krb_max.value();

dimensionedScalar fmmob = transportProperties.lookup("fmmob");
scalar fmmobValue = fmmob.value();

dimensionedScalar SF = transportProperties.lookup("SF");
scalar SFValue = SF.value();

dimensionedScalar sfbet = transportProperties.lookup("sfbet");
scalar sfbetValue = sfbet.value();

dimensionedScalar fmcap = transportProperties.lookup("fmcap");
scalar fmcapValue = fmcap.value();

dimensionedScalar epcap = transportProperties.lookup("epcap");
scalar epcapValue = epcap.value();

dimensionedScalar sigma_ba = transportProperties.lookup("sigma_ba");
scalar sigma_baValue = sigma_ba.value();

Info<< "Constants OK\n" << endl;

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

// volVectorField U
// (
//     IOobject
//     (
//         "U",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );

// #include "createPhi.H"

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

Info<< "Fluxes and velocities OK\n" << endl;

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
kra = correct_kra(Sb,Sa_minValue,Sb_minValue,kra_maxValue,a_expValue,fmmobValue,SFValue,sfbetValue,U,fmcapValue,epcapValue,mu_b,sigma_baValue);
surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));
Info<< "kra OK \n" << endl;
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
krb = correct_krb(Sb,Sa_minValue,Sb_minValue,krb_maxValue,b_expValue);
surfaceScalarField krbf ("krbf",fvc::interpolate(krb,"krb"));

//- mobility and fractional flow 
surfaceScalarField Maf ("Maf",Kf*kraf/mu_a);
surfaceScalarField Mbf ("Mbf",Kf*krbf/mu_b);
surfaceScalarField Mf ("Mf",Maf+Mbf);
surfaceScalarField Fbf ("Fbf",Mbf/Mf);

// volScalarField Fb ("Fb",(krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) ));
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

Info<< "Relative permeabilities OK\n" << endl;

Info<< "\nReading g" << endl;
dimensionedVector g_vector
(
    transportProperties.lookup("g_vector")
);

// volVectorField g
// (
//     IOobject
//     (
//         "g",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     g_vector
// );
// surfaceVectorField gf ("gf",fvc::interpolate(g,"g"));

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

//

// #include "readGravitationalAcceleration.H"

//

surfaceScalarField Laf ("Laf",rho_a*Kf*kraf/mu_a);
surfaceScalarField Lbf ("Lbf",rho_b*Kf*krbf/mu_b);
surfaceScalarField Lf ("Lf",Laf+Lbf);

surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());


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
Info<< "\nInitializing Pc" << endl;
pc = correct_pc(Sb);

Info<< "\nInitializing dpcdS" << endl;
volScalarField dpcdS
(
    IOobject
    (
        "dpcdS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);
Info<< "\nInitializing dpcdS" << endl;
dpcdS = correct_dpcdS(Sb);
Info<< "\nInitializing dpcdS" << endl;
surfaceScalarField dpcdSf ("dpcdSf",fvc::interpolate(dpcdS,"dpcdS"));

Info<< "\nInitializing phiPc" << endl;
surfaceScalarField phiPc("phiPc",Mbf * dpcdSf * fvc::snGrad(Sb) * mesh.magSf());


// SOURCE TERM

volScalarField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimless,0)
);

volScalarField qw_inj
(
    IOobject
    (
        "qw_inj",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    max(Q, scalar(0))  
);
qw_inj.write();

volScalarField qw_prod
(
    IOobject
    (
        "qw_prod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    min(Q, scalar(0))  
);
qw_prod.write();

std::remove("production.txt");
std::ofstream outFile("production.txt", std::ios::app);
outFile << "TimeStep,qw_inj,qw_prod,qo_prod\n"; 
outFile.close();

#include "createFvModels.H"
// #include "createFvConstraints.H"
