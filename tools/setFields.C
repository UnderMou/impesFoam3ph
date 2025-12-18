// TODO: Not ready yet! 
// Necessary to implement an OpenFOAM callable functionality that:
// > Read the mesh and return the cell index and x,y,z coordinates

// INITIALIZE K
forAll(K, i)
{   
    const scalar x = mesh.C()[i][0];
    const scalar y = mesh.C()[i][1];

    // K[i] = 1.0e-12 + 1.0e-13*x - 1.0e-13*y;
    // Info << i << "," << K[i] << "," << x << "," << y << endl;

    Info << i << "," << x << "," << y << endl;
}  