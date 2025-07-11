#include <cmath>

void printField(volScalarField field){
    forAll(field, i){
        Pout << field[i] << " ";
    }
    Pout << endl;
}

volScalarField correct_kra(const volScalarField Sa, 
                           const volScalarField Sb,
                           const volScalarField p,
                           const volScalarField K,
                           scalar Sa_min,
                           scalar Sb_min,
                           scalar Sc_min,
                           scalar kra_max,
                           scalar a,
                           scalar fmmob,
                           scalar SF,
                           scalar sfbet,
                           const volVectorField U,
                           scalar fmcap,
                           scalar epcap,
                           dimensionedScalar mu_b,
                           scalar sigma_ba,
                           scalar fmoil,
                           scalar floil,
                           scalar epoil)
{
    // // OPENFOAM ALIKE
    // volScalarField Se_a = (Sa-Sa_min)/(1.0-Sa_min-Sb_min-Sc_min);
    // volScalarField kra = kra_max * Foam::pow(Se_a, a);

    // FOSSIL ALIKE
    volScalarField kra = Sa;    // just for initialization
    scalar Se = 0.0;            // just for initialization
    scalar zero = 1e-8;
    forAll(kra,i){
        if (Sa[i]<=Sa_min + zero){
            kra[i]=0.0;
        }
        else if (Sa[i]>=1.0 - Sb_min - Sc_min - zero){
            kra[i]=kra_max;
        }
        else{
            Se = (Sa[i]-Sa_min)/(1.0-Sa_min-Sb_min-Sc_min);
            kra[i] = kra_max * Foam::pow(Se, a);
        }
    }

    // STARS Foam Model
    const double pi = 3.141592653589793; 

    // Fdry
    volScalarField Fdry = 0.5 + (1.0 / pi) * Foam::atan(sfbet * (Sb - SF));

    // // Fshear
    // volScalarField Fshear = Fdry;   // Just for initialization!
    // volScalarField Nca = (mu_b * mag(U)) / sigma_ba;
    // forAll(Nca,i){
    //     if (Nca[i] >= fmcap){
    //         Fshear[i] = Foam::pow(fmcap/Nca[i], epcap);
    //     }
    //     else{
    //         Fshear[i] = 1.0;
    //     }
    // } 

    // Fshear
    volScalarField Fshear = Fdry;   // Just for initialization!
    volVectorField gradP = fvc::grad(p);
    volScalarField Nca = (K * mag(gradP)) / sigma_ba;
    forAll(Nca,i){
        if (Nca[i] > fmcap){
            Fshear[i] = Foam::pow(fmcap/Nca[i], epcap);
        }
        else{
            Fshear[i] = 1.0;
        }
    }

    // // Fshear
    // volScalarField Fshear = Fdry;   // Just for initialization! 
    // forAll(Fshear,i){
    //     Fshear[i] = 1.0;
    // }

    // Foil
    volScalarField Foil = Fdry;   // Just for initialization!
    volScalarField Sc = 1.0 - Sa - Sb;
    forAll(Sc,i){
        if(fmoil < Sc[i] && Sc[i] <= 1.0 - Sa_min - Sb_min){
            Foil[i] = 0.0;
        }
        else if(floil <= Sc[i] && Sc[i] <= fmoil){
            Foil[i] = Foam::pow((fmoil-Sc[i])/(fmoil-floil),epoil);
        }
        else{
            Foil[i] = 1.0;
        }
    }

    // MRF
    volScalarField MRF = 1.0 / (1.0 + fmmob*Fdry*Fshear*Foil);
    // volScalarField MRF = 1.0 / (1.0 + fmmob*Fdry*Foil);
    kra *= MRF;

    return kra;
}

volScalarField correct_krb(const volScalarField Sb,
                           scalar Sa_min,
                           scalar Sb_min,
                           scalar Sc_min,
                           scalar krb_max,
                           scalar b)
{
    // OPENFOAM ALIKE
    // volScalarField Se_b = (Sb-Sb_min)/(1.0-Sa_min-Sb_min-Sc_min);
    // volScalarField krb = krb_max * Foam::pow(Se_b, b); 

    // FOSSIL ALIKE
    volScalarField krb = Sb;    // just for initialization
    scalar Se = 0.0;            // just for initialization
    scalar zero = 1e-8;
    forAll(krb,i){
        if (Sb[i]<=Sb_min + zero){
            krb[i]=0.0;
        }
        else if (Sb[i]>=1.0 - Sa_min - Sc_min - zero){
            krb[i]=krb_max;
        }
        else{
            Se = (Sb[i]-Sb_min)/(1.0-Sa_min-Sb_min-Sc_min);
            krb[i] = krb_max * Foam::pow(Se, b);
        }
    }

    return krb;
}

volScalarField correct_krc(const volScalarField Sa,
                           const volScalarField Sb,
                           scalar Sa_min,
                           scalar Sb_min,
                           scalar Sc_min,
                           scalar krc_max,
                           scalar c)
{
    // OPENFOAM ALIKE
    // volScalarField Se_a = (Sa-Sa_min)/(1.0-Sa_min-Sb_min-Sc_min);
    // volScalarField Se_b = (Sb-Sb_min)/(1.0-Sa_min-Sb_min-Sc_min);
    // volScalarField krc = krc_max * Foam::pow(scalar(1.0)-Se_a-Se_b, c); 

    // FOSSIL ALIKE
    volScalarField Sc = scalar(1.0) - Sa - Sb;
    volScalarField krc = Sc;    // just for initialization
    scalar Se = 0.0;            // just for initialization
    scalar zero = 1e-8;
    forAll(krc,i){
        if (Sc[i]<=Sc_min + zero){
            krc[i]=0.0;
        }
        else if (Sc[i]>=1.0 - Sb_min - Sa_min - zero){
            krc[i]=krc_max;
        }
        else{
            Se = (Sc[i]-Sc_min)/(1.0-Sa_min-Sb_min-Sc_min);
            krc[i] = krc_max * Foam::pow(Se, c);
        }
    }

    // volScalarField Sc = 1.0 - Sa - Sb;
    // volScalarField Se_c = (Sc-Sc_min)/(1.0-Sa_min-Sb_min-Sc_min);
    // volScalarField krc = krc_max * Foam::pow(Se_c, c); 

    return krc;
}

volScalarField correct_Nca(const volScalarField p,
                           const volScalarField K,
                           scalar sigma_ba)
{
    volVectorField gradP = fvc::grad(p);
    volScalarField Nca = (K * mag(gradP)) / sigma_ba;
    return Nca;
}

volScalarField correct_Fshear(const volScalarField p,
                              const volScalarField K,
                              scalar fmcap,
                              scalar epcap,
                              scalar sigma_ba)
{
    // Fshear
    volScalarField Fshear = K;   // Just for initialization!
    volVectorField gradP = fvc::grad(p);
    volScalarField Nca = (K * mag(gradP)) / sigma_ba;
    forAll(Nca,i){
        if (Nca[i] > fmcap){
            Fshear[i] = Foam::pow(fmcap/Nca[i], epcap);
        }
        else{
            Fshear[i] = 1.0;
        }
    }

    return Fshear;
}

volScalarField correct_Foil(const volScalarField Sa, 
                            const volScalarField Sb,
                            scalar Sa_min,
                            scalar Sb_min,
                            scalar fmoil,
                            scalar floil,
                            scalar epoil)
{
    // Foil
    volScalarField Foil = Sa;   // Just for initialization!
    volScalarField Sc = 1.0 - Sa - Sb;
    forAll(Sc,i){
        if(fmoil < Sc[i] && Sc[i] <= 1.0 - Sa_min - Sb_min){
            Foil[i] = 0.0;
        }
        else if(floil <= Sc[i] && Sc[i] <= fmoil){
            Foil[i] = Foam::pow((fmoil-Sc[i])/(fmoil-floil),epoil);
        }
        else{
            Foil[i] = 1.0;
        }
    }
    return Foil;
}