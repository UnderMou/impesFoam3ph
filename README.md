# impesFoam3ph
Implementation of the IMPES method on the 3-phase flow in porous media with foam.

This solver was based on the work: 
[An open-source toolbox for multiphase flow in porous media](https://www.sciencedirect.com/science/article/pii/S0010465514003403) P Horgue, C Soulaine, J Franc, R Guibert, G Debenest Computer Physics Communications 187, 217-226, 2015

New features added:
- Three-phase flow modeling
- Foam modeling:
   - Implicit-texture foam model added [(STARS foam model)](https://www.sciencedirect.com/science/article/pii/S1875510018300878).
   - Mechanistic foam model added [(Ashoori foam model)](https://www.sciencedirect.com/science/article/pii/S0927775711000124).
- Surfactant concentration transport PDE.
  
---
**IMPORTANT**: This repository is part of ongoing research associated with a scientific paper that is currently under preparation/submission. The contents are under active development and may change.

![Under Development](https://img.shields.io/badge/status-under%20development-orange)
___

## Instructions

### 1. Installing OpenFOAM on Your Machine
To use the *impesFoam3ph* solver, you need to install OpenFOAM.org v9 on your machine (Linux).

For Ubuntu systems, installation can be done by following this tutorial: [Ubuntu Installation Guide](https://openfoam.org/download/9-ubuntu/).

With OpenFOAM properly installed on your machine according to the tutorial above, you can proceed to compile the *impesFoam3ph* solver.

### 2. Compiling the Files

   2.1 **Compiling the `impesFoam3ph` Solver**  

   For now, there are two solvers implemented:
   - impesFoam3ph for the 3-phase cases
   - impesFoam2ph for the 2-phase cases

   Navigate to your `$FOAM_RUN` directory, and clone this repository there by:
   ```bash
   cd $FOAM_RUN
   git clone <repository-link>
   ```

   Once the repository is cloned, go to the repository directory and compile the solvers by executing:
   ```bash
   cd <the-directory-of-the-cloned-repository>'
   chmod +x build_solvers.sh
   ./build_solvers.sh
   ```

   Verify that the solver compiled successfully by listing the contents of the `$FOAM_USER_APPBIN` directory:
   ```bash
   ls $FOAM_USER_APPBIN
   ```
   If `impesFoam3ph` and `impesFoam2ph` appear in the listing, the compilation was successful.

   2.2 **Compiling the `myPorousMediumBCs` Library**  
   The `myPorousMediumBCs` directory contains custom boundary conditions for porous media applications (for now, there is just `gradPressureDarcy`). This library must also be compiled.

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

   . for 3-phase cases
   ```bash
   blockMesh
   impesFoam3ph
   ```
   . for 2-phase cases
   ```bash
   blockMesh
   impesFoam2ph
   ```

   To view the results, open ParaView and load the `field.foam` file located in the case directory.

