volScalarField correct_Sb(volScalarField Sb, scalar Sa_min, scalar Sb_min){
    // double Sb_min = 0.15;
    // double Sa_min = 0.1;

    // forAll(Sb, i){
    //     if (Sb[i] <= Sb_min) Sb[i] = Sb_min;
    //     else if (Sb[i] >= 1.0-Sa_min) Sb[i] = 1.0-Sa_min;
    // }

    return Sb;
}