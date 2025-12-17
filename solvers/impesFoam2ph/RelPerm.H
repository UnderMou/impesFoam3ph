void correct_kra
(
    volScalarField& kra,
    const volScalarField Sb, 
    const scalar Sa_min,
    const scalar Sb_min,
    const scalar kra_max,
    const scalar a
)
{
    volScalarField Se = (Sb-Sb_min)/(1.0-Sa_min-Sb_min);
    kra = kra_max * Foam::pow(scalar(1.0)-Se, a);
}

void correct_krb
(
    volScalarField& krb,
    const volScalarField Sb,
    const scalar Sa_min,
    const scalar Sb_min,
    const scalar krb_max,
    const scalar b
)
{
    volScalarField Se = (Sb-Sb_min)/(1.0-Sa_min-Sb_min);
    krb = krb_max * Foam::pow(Se, b); 
}


