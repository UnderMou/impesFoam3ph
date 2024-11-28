# impesFoam3ph
Implementation of the IMPES method on the 3 phase flow in porous media with foam.

## Instructions

### 1. Installing OpenFOAM on Your Machine
To use the *impesFoam3ph* solver, you need to install OpenFOAM.org v9 on your machine (Linux).

For Ubuntu systems, installation can be done by following this tutorial: [Ubuntu Installation Guide](https://openfoam.org/download/9-ubuntu/).

With OpenFOAM properly installed on your machine according to the tutorial above, you can proceed to compile the *impesFoam3ph* solver.

### 2. Compiling the Files

   2.1 **Compiling the `impesFoam3ph` Solver**  
   Navigate to your `$FOAM_RUN` directory, and clone this repository there by:
   ```bash
   cd $FOAM_RUN
   git clone <repository-link>
   ```

   Once the repository is cloned, go to the repository directory and compile the `impesFoam3ph` solver by executing:
   ```bash
   cd <the-directory-of-the-cloned-repository>'
   wmake
   ```

   Verify that the solver compiled successfully by listing the contents of the `$FOAM_USER_APPBIN` directory:
   ```bash
   ls $FOAM_USER_APPBIN
   ```
   If `impesFoam3ph` appears in the listing, the compilation was successful.

   2.2 **Compiling the `myPorousMediumBCs` Library**  
   The `myPorousMediumBCs` directory contains custom boundary conditions for porous media applications (for now there is just `gradPressureDarcy`). This library must also be compiled.

   Navigate to the `myPorousMediumBCs` directory, and compile it by running:
   ```bash
   wmake libso .
   ```

   Check if the library was successfully compiled by listing the contents of the `$FOAM_USER_LIBBIN` directory:
   ```bash
   ls $FOAM_USER_LIBBIN
   ```
   If `libmyPorousMediumBCs.so` is listed, the compilation was successful.

### 3. Running an Example

   To verify that everything is set up correctly, go to any case within the *Examples* directory and run it using the following commands:
   ```bash
   blockMesh
   impesFoam3ph_v1
   ```
   `blockMesh` is the OpenFOAM meshing tool, and `impesFoam3ph_v1` will initialize and run the solver.

   To view the results, open ParaView and load the `field.foam` file located in the case directory.

