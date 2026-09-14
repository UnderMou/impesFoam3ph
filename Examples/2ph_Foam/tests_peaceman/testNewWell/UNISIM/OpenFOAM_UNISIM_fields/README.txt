UNISIM OpenFOAM mapped fields

Mesh cells: 38208

Files:
  PORO  : original PORO values, dimensionless.
  PERMX : original PERMX values converted from mD to m^2.
  K     : same SI permeability values as PERMX, provided because the OpenFOAM solver uses field name K.

Conversion:
  1 mD = 9.869233000e-16 m^2

Statistics:
  PORO min/max   = 0 / 0.3
  PERMX mD min/max = 1 / 1190
  PERMX m^2 min/max = 9.86923300e-16 / 1.17443873e-12

Boundary patch:
  defaultFaces -> zeroGradient

Usage:
  Copy PORO and K (or PERMX) into the case 0/ directory after converting
  UNISIM-I-D_clean_compact.vtk with vtkUnstructuredToFoam.

The mapping follows the cell ordering in UNISIM-I-D_clean_compact.vtk.
