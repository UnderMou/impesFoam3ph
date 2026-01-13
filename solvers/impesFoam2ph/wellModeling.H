void correct_qb
(
    volScalarField& qb,
    const volScalarField& Fb, 
    const volScalarField& qt_inj,
    const volScalarField& qt_prod
)   
{   
    const scalar Fb_inj = 1.0;

    qb = Fb_inj * qt_inj - Fb * qt_prod;

}