#include <cmath>

void printField(volScalarField field){
    forAll(field, i){
        Pout << field[i] << " ";
    }
    Pout << endl;
}

volScalarField correct_kra(const volScalarField Sa, 
                           const volScalarField Sb,
                           const volScalarField Sc,
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
    // scalar zero = 1e-8;
    // volScalarField Se_a = (Sa-Sa_min+zero)/(1.0-Sa_min-Sb_min-Sc_min);
    // volScalarField kra = kra_max * Foam::pow(Se_a, a);

    volScalarField kra = Sa;    // just for initialization
    scalar Se = 0.0;            // just for initialization
    scalar zero = 1e-8;
    scalar Sa_max = scalar(1.0) - Sb_min - Sc_min;
    forAll(kra,i){
        if (Sa[i]<=Sa_min + zero){
            kra[i]=0.0;
        }
        else if (Sa[i]>=1.0 - (Sa_max + zero)){
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

    // // Fshear
    // volScalarField Fshear = Fdry;   // Just for initialization!
    // volVectorField gradP = fvc::grad(p);
    // volScalarField Nca = (K * mag(gradP)) / sigma_ba;
    
    // forAll(Nca,i){
    //     if (Nca[i] > fmcap){
    //         Fshear[i] = Foam::pow(fmcap/max(Nca[i],zero), epcap);
    //     }
    //     else{
    //         Fshear[i] = 1.0;
    //     }
    // }

    // // Fshear
    // volScalarField Fshear = Fdry;   // Just for initialization! 
    // forAll(Fshear,i){
    //     Fshear[i] = 1.0;
    // }

    // Foil
    volScalarField Foil = Fdry;   // Just for initialization!
    // volScalarField Sc = 1.0 - Sa - Sb;
    // forAll(Sc,i){
    //     if(fmoil < Sc[i] && Sc[i] <= 1.0 - Sa_min - Sb_min + zero){
    //         Foil[i] = 0.0;
    //     }
    //     else if(floil <= Sc[i] && Sc[i] <= fmoil){
    //         Foil[i] = Foam::pow((fmoil-Sc[i])/(fmoil-floil),epoil);
    //     }
    //     else{
    //         Foil[i] = 1.0;
    //     }
    // }
    forAll(Sc,i){
        if(Sc[i] <= floil + zero){
            Foil[i] = 1.0;
        }
        else if(Sc[i] >= fmoil - zero){
            Foil[i] = 0.0;
        }
        else{
            Foil[i] = Foam::pow((fmoil-Sc[i])/(fmoil-floil),epoil);
        }
    }

    // MRF
    // volScalarField MRF = 1.0 / (1.0 + max(fmmob*Fdry*Fshear*Foil,zero));
    volScalarField MRF = 1.0 / (1.0 + fmmob*Fdry*Foil);
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
    // scalar zero = 1e-8;
    // volScalarField Se_b = (Sb-Sb_min+zero)/(1.0-Sa_min-Sb_min-Sc_min+VSMALL);
    // volScalarField krb = krb_max * Foam::pow(Se_b, b); 

    volScalarField krb = Sb;    // just for initialization
    scalar Se = 0.0;            // just for initialization
    scalar zero = 1e-8;
    scalar Sb_max = scalar(1.0) - Sa_min - Sc_min;
    forAll(krb,i){
        if (Sb[i]<=Sb_min + zero){
            krb[i]=0.0;
        }
        else if (Sb[i]>=1.0 - (Sb_max + zero)){
            krb[i]=krb_max;
        }
        else{
            Se = (Sb[i]-Sb_min)/(1.0-Sa_min-Sb_min-Sc_min);
            krb[i] = krb_max * Foam::pow(Se, b);
        }
    }


    return krb;
}

volScalarField correct_krc(const volScalarField Sc,
                           scalar Sa_min,
                           scalar Sb_min,
                           scalar Sc_min,
                           scalar krc_max,
                           scalar c)
{

    // scalar zero = 1e-8;
    // volScalarField Se_c = (Sc-Sc_min)/(1.0-Sa_min-Sb_min-Sc_min);
    // volScalarField krc = krc_max * Foam::pow(Se_c, c); 

    volScalarField krc = Sc;    // just for initialization
    scalar Se = 0.0;            // just for initialization
    scalar zero = 1e-8;
    scalar Sc_max = scalar(1.0) - Sa_min - Sb_min;
    forAll(krc,i){
        if (Sc[i]<=Sc_min + zero){
            krc[i]=0.0;
        }
        else if (Sc[i]>=1.0 - (Sc_max + zero)){
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

volScalarField correct_Sb(volScalarField Sb, scalar Sa_min, scalar Sb_min, scalar Sc_min){
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    scalar zero = 1e-4;

    // forAll(Sb, i){
    //     if (Sb[i] <= (Sb_min + zero)){
    //         Sb[i] = Sb_min + zero;
    //     }
    //     else if (Sb[i] >= (1.0 - Sa_min - Sc_min - zero)){
    //         Sb[i] = 1.0 - Sa_min - Sc_min - zero;
    //     }
    // }

    forAll(Sb, i){
        if (Sb[i] <= (zero)){
            Sb[i] = zero;
        }
        else if (Sb[i] >= (1.0 - zero)){
            Sb[i] = 1.0 - zero;
        }
    }

    return Sb;
}

volScalarField correct_Sa(volScalarField Sa, scalar Sa_min, scalar Sb_min, scalar Sc_min){
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    scalar zero = 1e-4;

    // forAll(Sa, i){
    //     if (Sa[i] <= (Sa_min + zero)){
    //         Sa[i] = Sa_min + zero;
    //     }
    //     else if (Sa[i] >= (1.0 - Sb_min - Sc_min - zero)){
    //         Sa[i] = 1.0 - Sb_min - Sc_min - zero;
    //     }
    // }

    forAll(Sa, i){
        if (Sa[i] <= (zero)){
            Sa[i] = zero;
        }
        else if (Sa[i] >= (1.0 - zero)){
            Sa[i] = 1.0 - zero;
        }
    }

    return Sa;
}

volScalarField correct_Sc(volScalarField Sc, scalar Sa_min, scalar Sb_min, scalar Sc_min){
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    scalar zero = 1e-4;

    // forAll(Sc, i){
    //     if (Sc[i] <= (Sc_min + zero)){
    //         Sc[i] = Sc_min + zero;
    //     }
    //     else if (Sc[i] >= (1.0 - Sb_min - Sa_min - zero)){
    //         Sc[i] = 1.0 - Sb_min - Sa_min - zero;
    //     }
    // }

    forAll(Sc, i){
        if (Sc[i] <= (zero)){
            Sc[i] = zero;
        }
        else if (Sc[i] >= (1.0 - zero)){
            Sc[i] = 1.0 - zero;
        }
    }

    return Sc;
}