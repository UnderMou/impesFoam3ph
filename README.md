# impesFoam3ph
Computational implementation of the IMplicit-Pressure Explicit-Saturation (IMPES) method on the three-phase flow in porous media with foam using OpenFOAM framework.

Main features:
- Three-phase flow modeling
- Foam modeling:
   - Implicit-texture foam model added [(STARS foam model)](https://www.sciencedirect.com/science/article/pii/S1875510018300878).
   - Mechanistic foam model added [(Ashoori foam model)](https://www.sciencedirect.com/science/article/pii/S0927775711000124).
- Surfactant concentration transport in aqueous phase PDE. [(Surfactant Concentration transport in aqueous phase)](https://epubs.siam.org/doi/10.1137/23M1566649) accounting with surfactant adsorption on solid phase.
  
---
**IMPORTANT**: This repository is part of an ongoing research associated with a scientific paper that is currently under preparation/submission. The contents are under active development and may change.

![Under Development](https://img.shields.io/badge/status-under%20development-orange)
___

## Examples

This repository provides a collection of example cases covering a wide range of multiphase flow scenarios in porous media, including capillary effects, gravity, foam, and well-driven flows.

### Folder Structure

### Two-Phase Flow

- **2ph_CapillaryPressure**
  - `capGrav_equi`: Vertical gas–water column reaching capillary–gravity equilibrium.

- **2ph_Foam**
  - `2Layers_long`: Co-injection in horizontally layered heterogeneous medium.
  - `2Layers_Seq`: Sequential (flow-aligned) heterogeneity with vertical zones.
  - `ashoori_diff_Kcs`: Sensitivity study on foam generation rate ($K_c$).
  - `ashoori_paper`: Base case inspired by Ashoori et al. foam model.
  - `ChessField`: Checkerboard permeability field with strong heterogeneity.
  - `Core16_foamQualityScan`: Foam quality scan with experimental parameter fitting.
  - `Core21`: Core-scale foam simulation case.
  - `SPE10_L36`: Full heterogeneous field from SPE10 Layer 36.
  - `SurfactantTransport`: Foam with surfactant transport and adsorption effects.

- **2ph_Gravity**
  - `2D_homogeneous`: Vertical homogeneous domain with gravity and weak foam.
  - `Durlofsky_Abreu`: 1D counter-gravity injection benchmark.

- **2ph_Well**
  - `5spot_panday`: Five-spot pattern water injection for oil displacement.

---

### Three-Phase Flow

- **3ph_CapillaryPressure**
  - `LeveretLewis_EduardoAbreu`: Riemann problems with capillary diffusion (Leverett–Lewis model).

- **3ph_Foam**
  - `LozanoFoam`: Three-phase foam model based on Riemann problem formulation.
  - `LyuTang_CTscan`: Foam simulation inspired by CT-scan experiments.
  - `LyuTang_fracFlow_verification`: Fractional flow validation case for the three-phase flow.
  - `Mehrabi_surfTransp`: Secondary and tertiary recovery with surfactant transport.

- **3ph_Gravity**
  - `EduardoAbreu`: 1D counter-gravity injection three-phase flow benchmark.

- **3ph_Well**
  - `3D_homogeneous`: Three-phase flow in homogeneous reservoir with wells.
  - `3D_SPE10field_10layers`: Multilayer SPE10-based heterogeneous reservoir simulation.


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

   Once the repository is cloned, go to the repository directory and compile the files by executing:
   ```bash
   cd <the-directory-of-the-cloned-repository>'
   chmod +x compile.sh
   ./compile.sh
   ```

   2.2 **Checking the compiled files** 

   Verify that the solver compiled successfully by listing the contents of the `$FOAM_USER_APPBIN` directory:
   ```bash
   ls $FOAM_USER_APPBIN
   ```
   If `impesFoam3ph` and `impesFoam2ph` appear in the listing, the compilation was successful.

   Check if the libraries was successfully compiled by listing the contents of the `$FOAM_USER_LIBBIN` directory:
   ```bash
   ls $FOAM_USER_LIBBIN
   ```
   If all the following files are listed, the compilation was successfull:
   `libboundaryConditions_2ph.so`, `libfoamModels_2ph.so`, `librelativePermeabilityModels_2ph.so`, `libwellModels_2ph.so`, `libboundaryConditions_3ph.so`, `libfoamModels_3ph.so`, `librelativePermeabilityModels_3ph.so`, `libwellModels_3ph.so`, `libcapillaryPressureModels_2ph.so`, `libisothermModels.so` and `libsurfactantTransportModels.s`.

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

## References
This solver was based on the work: 

[An open-source toolbox for multiphase flow in porous media](https://www.sciencedirect.com/science/article/pii/S0010465514003403) P Horgue, C Soulaine, J Franc, R Guibert, G Debenest Computer Physics Communications 187, 217-226, 2015
