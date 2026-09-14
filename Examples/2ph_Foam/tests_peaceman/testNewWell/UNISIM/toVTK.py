import meshio

mesh = meshio.read("UNISIM-I-D_optionA.vtu")

meshio.write(
    "UNISIM-I-D_optionA.vtk",
    mesh,
    file_format="vtk42",
    binary=False
)