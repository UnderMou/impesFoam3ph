void printField(volScalarField field){
    forAll(field, i){
        Pout << field[i] << " ";
    }
    Pout << endl;
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