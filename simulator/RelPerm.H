#include <cmath>

void printField(volScalarField field){
    forAll(field, i){
        Pout << field[i] << " ";
    }
    Pout << endl;
}

volScalarField correct_kra(const volScalarField Sa, 
                           const volScalarField Sb,
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

    // Foil
    volScalarField Sc = 1.0 - Sa - Sb;
    // volScalarField Foil = Fdry;   // Just for initialization!
    // forAll(Sc,i){
    //     if(fmoil <= Sc[i] && Sc[i] <= 1.0 - Sa_min - Sb_min){
    //         Foil[i] = 0.0;
    //     }
    //     else if(floil <= Sc[i] && Sc[i] < fmoil){
    //         Foil[i] = Foam::pow((fmoil-Sc[i])/(fmoil-floil),epoil);
    //     }
    //     else{
    //         Foil[i] = 1.0;
    //     }
    // }

    volScalarField Foil = 0.5 + (1.0 / pi) * Foam::atan(epoil * (fmoil - Sc));

    // MRF
    // volScalarField MRF = 1.0 / (1.0 + fmmob*Fdry*Fshear);
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

    volScalarField Se_b = (Sb-Sb_min)/(1.0-Sa_min-Sb_min-Sc_min);

    volScalarField krb = krb_max * Foam::pow(Se_b, b); 

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

    return krc;
}