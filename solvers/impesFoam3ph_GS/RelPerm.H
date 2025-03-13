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
    volScalarField Se_a = (Sa-Sa_min)/(1.0-Sa_min-Sb_min-Sc_min);

    volScalarField kra = kra_max * Foam::pow(Se_a, a);

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

    volScalarField Se_b = (Sb-Sb_min)/(1.0-Sa_min-Sb_min-Sc_min);
    volScalarField krb = krb_max * Foam::pow(Se_b, b); 

    // 2-spot case from Carolina Xavier
    // volScalarField krb = Se_b;

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

    volScalarField Se_a = (Sa-Sa_min)/(1.0-Sa_min-Sb_min-Sc_min);
    volScalarField Se_b = (Sb-Sb_min)/(1.0-Sa_min-Sb_min-Sc_min);
    volScalarField krc = krc_max * Foam::pow(scalar(1.0)-Se_a-Se_b, c); 

    // 2-spot case from Carolina Xavier
    // volScalarField krc = krc_max * Foam::pow(scalar(1.0)-Se_b, c); 


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

volScalarField correct_Sb(volScalarField Sb, scalar Sa_min, scalar Sb_min, scalar Sc_min){
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    scalar zero = 1e-8;

    forAll(Sb, i){
        if (Sb[i] <= (Sb_min + zero)){
            Sb[i] = Sb_min + zero;
        }
        else if (Sb[i] >= (1.0 - Sa_min - Sc_min - zero)){
            Sb[i] = 1.0 - Sa_min - Sc_min - zero;
        }
    }

    return Sb;
}

volScalarField correct_Sa(volScalarField Sa, scalar Sa_min, scalar Sb_min, scalar Sc_min){
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    scalar zero = 1e-8;

    forAll(Sa, i){
        if (Sa[i] <= (Sa_min + zero)){
            Sa[i] = Sa_min + zero;
        }
        else if (Sa[i] >= (1.0 - Sb_min - Sc_min - zero)){
            Sa[i] = 1.0 - Sb_min - Sc_min - zero;
        }
    }

    return Sa;
}