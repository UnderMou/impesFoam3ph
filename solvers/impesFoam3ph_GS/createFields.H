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

dimensionedScalar Sa_min = transportProperties.lookup("Sa_min");
scalar Sa_minValue = Sa_min.value();

dimensionedScalar Sb_min = transportProperties.lookup("Sb_min");
scalar Sb_minValue = Sb_min.value();

dimensionedScalar Sc_min = transportProperties.lookup("Sc_min");
scalar Sc_minValue = Sc_min.value();

dimensionedScalar a_exp = transportProperties.lookup("a_exp");
scalar a_expValue = a_exp.value();

dimensionedScalar b_exp = transportProperties.lookup("b_exp");
scalar b_expValue = b_exp.value();

dimensionedScalar c_exp = transportProperties.lookup("c_exp");
scalar c_expValue = c_exp.value();

dimensionedScalar kra_max = transportProperties.lookup("kra_max");
scalar kra_maxValue = kra_max.value();

dimensionedScalar krb_max = transportProperties.lookup("krb_max");
scalar krb_maxValue = krb_max.value();

dimensionedScalar krc_max = transportProperties.lookup("krc_max");
scalar krc_maxValue = krc_max.value();

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

dimensionedScalar fmoil = transportProperties.lookup("fmoil");
scalar fmoilValue = fmoil.value();

dimensionedScalar floil = transportProperties.lookup("floil");
scalar floilValue = floil.value();

dimensionedScalar epoil = transportProperties.lookup("epoil");
scalar epoilValue = epoil.value();

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
kra = correct_kra(Sa,Sb,p,K,Sa_minValue,Sb_minValue,Sc_minValue,kra_maxValue,a_expValue,fmmobValue,SFValue,sfbetValue,U,fmcapValue,epcapValue,mu_b,sigma_baValue,fmoilValue,floilValue,epoilValue);
surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));

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
krb = correct_krb(Sb,Sa_minValue,Sb_minValue,Sc_minValue,krb_maxValue,b_expValue);
surfaceScalarField krbf ("krbf",fvc::interpolate(krb,"krb"));

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
krc = correct_krc(Sa,Sb,Sa_minValue,Sb_minValue,Sc_minValue,krc_maxValue,c_expValue);
surfaceScalarField krcf ("krcf",fvc::interpolate(krc,"krc"));

//- mobility and fractional flow 
surfaceScalarField Maf ("Maf",Kf*kraf/mu_a);
surfaceScalarField Mbf ("Mbf",Kf*krbf/mu_b);
surfaceScalarField Mcf ("Mcf",Kf*krcf/mu_c);
surfaceScalarField Mf ("Mf",Maf+Mbf+Mcf);
surfaceScalarField Faf ("Faf",Maf/Mf);
surfaceScalarField Fbf ("Fbf",Mbf/Mf);
surfaceScalarField Fcf ("Fcf",Mcf/Mf);

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
Fb = (krb/mu_b) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );

volScalarField Fc
(
    IOobject
    (
        "Fc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);
Fc = (krc/mu_c) / ( (kra/mu_a) + (krb/mu_b) + (krc/mu_c) );

// GRAVITY EFFECTS

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

surfaceScalarField Lacf ("Lacf",Fcf*Laf - Faf*Lcf);
surfaceScalarField Labf ("Labf",Fbf*Laf - Faf*Lbf);
surfaceScalarField Lbcf ("Lbcf",Fcf*Lbf - Fbf*Lcf);
surfaceScalarField Lbaf ("Lbaf",Faf*Lbf - Fbf*Laf);

surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());

////////////////////////////////////////////////////////

volScalarField Nca
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
    0.0
);
Nca = correct_Nca(p, K, sigma_baValue);

volScalarField Fshear
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
    0.0
);
Fshear = correct_Fshear(p, K, fmcapValue, epcapValue, sigma_baValue);

volScalarField Foil
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
    0.0
);
Foil = correct_Foil(Sa, Sb, Sa_minValue, Sb_minValue, fmoilValue, floilValue, epoilValue);

////////////////////////////////////////////////////////



std::remove("production.txt");
std::ofstream outFile("production.txt", std::ios::app);
outFile << "TimeStep,CellID,qw,qo\n"; 
outFile.close();

#include "createFvModels.H"
// #include "createFvConstraints.H"