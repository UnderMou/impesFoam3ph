#include <cmath>

void printField(volScalarField field){
    forAll(field, i){
        Pout << field[i] << " ";
    }
    Pout << endl;
}

volScalarField correct_kra(const volScalarField Sb, 
                           scalar Sa_min,
                           scalar Sb_min,
                           scalar kra_max,
                           scalar a,
                           scalar fmmob,
                           scalar SF,
                           scalar sfbet,
                           const volVectorField U,
                           scalar fmcap,
                           scalar epcap,
                           dimensionedScalar mu_b,
                           scalar sigma_ba)
{
    // parameters 
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    // double kra_max = 0.5;
    // double a = 8.0;

    volScalarField Se = (Sb-Sb_min)/(1.0-Sa_min-Sb_min);

    volScalarField kra = kra_max * Foam::pow(scalar(1.0)-Se, a);

    // forAll(kra,i){
    //     if (Sb[i]<=Sb_min){
    //         kra[i]=kra_max;
    //     }
    //     else if (Sb[i]>=1.0-Sb_min){
    //         kra[i]=0.0;
    //     }
    //     else{
    //         if (kra[i]>=kra_max) kra[i]=kra_max;
    //     }
    // }

    const double pi = 3.141592653589793; 

    // double fmmob = 170.113;
    // double SF = 0.45656;
    // double sfbet = 504.143;

    volScalarField Fdry = 0.5 + (1.0 / pi) * Foam::atan(sfbet * (Sb - SF));

    volScalarField Fshear = Fdry;   // Just for initialization!
    volScalarField Nca = (mu_b * mag(U)) / sigma_ba;
    forAll(Nca,i){
        if (Nca[i] >= fmcap){
            Fshear[i] = Foam::pow(fmcap/Nca[i], epcap);
        }
        else{
            Fshear[i] = 1.0;
        }
    } 

    volScalarField MRF = 1.0 / (1.0 + fmmob*Fdry*Fshear);
    kra *= MRF;

    return kra;
}

volScalarField correct_krb(const volScalarField Sb, scalar Sa_min, scalar Sb_min, scalar krb_max, scalar b){
    // parameters 
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    // double krb_max = 1.0;
    // double b = 3.0;

    volScalarField Se = (Sb-Sb_min)/(1.0-Sa_min-Sb_min);

    volScalarField krb = krb_max * Foam::pow(Se, b); 

    // forAll(krb,i){
    //     if (Sb[i]<=Sb_min){
    //         krb[i]=0.0;
    //     }
    //     else if (Sb[i]>=1.0-Sa_min){
    //         krb[i]=krb_max;
    //     }
    //     else{
    //         if (krb[i]>=krb_max) krb[i]=krb_max;
    //     }
    // }

    return krb;
}


volScalarField correct_Sb(volScalarField Sb, scalar Sa_min, scalar Sb_min){
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;

    // forAll(Sb, i){
    //     if (Sb[i] <= Sb_min) Sb[i] = Sb_min;
    //     else if (Sb[i] >= 1.0-Sa_min) Sb[i] = 1.0-Sa_min;
    // }

    return Sb;
}