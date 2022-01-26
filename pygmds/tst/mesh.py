from pygmds import *

mod = MeshModel(DIM3|F|E|N|F2N|F2E|F2N|E2N|N2F)
m = Mesh(mod)
ios = IOService(m)

# File reading
reader = VTKReader(ios)
reader.setCellOptions(MeshModel(N|F))
reader.setDataOptions(MeshModel(N|F))
reader.read("quarter_crown.vtk")

# Mesh statistics
print("Nb nodes ", m.getNbNodes())
print("Nb faces ", m.getNbFaces())