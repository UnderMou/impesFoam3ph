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
                           scalar sigma_ba,
                           scalar fmsurf,
                           scalar epsurf,
                           const volScalarField Cs)
{
    // parameters 
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;
    // double kra_max = 0.5;
    // double a = 8.0;

    volScalarField Se = (Sb-Sb_min)/(1.0-Sa_min-Sb_min);
    // volScalarField Sa = 1.0 - Sb; 
    // volScalarField Se = (Sa-Sa_min)/(1.0-Sa_min-Sb_min);
    // std::cout.precision(8);
    

    // forAll(Se, i)
    // {   
    //     if (Se[i]>=1.0)
    //     {
    //         Se[i] = 1.0 - 1e-6;
    //     }
    //     else if (Se[i]<=0.0)
    //     {
    //         Se[i] = 1e-6;
    //     }
    //     // Info << Se[i] << " ";
    // }
    // Info << endl;

    // Info << "Se: " << " Min(Se) = " << gMin(Se) << " Max(Se) = " << gMax(Se) << endl;
     

    volScalarField kra = kra_max * Foam::pow(scalar(1.0)-Se, a);
    // volScalarField kra = kra_max * Foam::pow(Se, a);

    // Info << "ok" << endl;

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



    volScalarField Fsurf = Fdry;   // Just for initialization!
    forAll(Cs,i){
        if (Cs[i] >= fmsurf){
            Fsurf[i] = 1.0;
        }
        else{
            Fsurf[i] = Foam::pow(Cs[i]/fmsurf, epsurf);
        }
    }


    // volScalarField MRF = 1.0 / (1.0 + fmmob*Fdry);
    volScalarField MRF = 1.0 / (1.0 + fmmob*Fdry*Fshear*Fsurf);

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
    scalar zero = 1e-8;

    forAll(Sb, i){
        if (Sb[i] <= (Sb_min + zero)){
            Sb[i] = Sb_min + zero;
        }
        else if (Sb[i] >= (1.0 - Sa_min - zero)){
            Sb[i] = 1.0-Sa_min - zero;
        }
    }

    return Sb;
}

