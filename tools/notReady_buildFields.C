// This is not even close to be ready!
// Just an idea for initialization of fields. 

// INITIALIZE K
    forAll(K, i)
    {   
        const scalar x = mesh.C()[i][0];
        const scalar y = mesh.C()[i][1];

        // K[i] = 1.0e-12 + 1.0e-13*x - 1.0e-13*y;
        // Info << i << "," << K[i] << "," << x << "," << y << endl;

        Info << i << "," << x << "," << y << endl;
    }  