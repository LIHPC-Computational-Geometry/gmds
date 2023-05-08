#!/usr/bin/env python
import sys

import pygmds

file_name = str(sys.argv[1]) + "/quarter_crown.vtk"

print("file  :" + file_name)
mod = pygmds.MeshModel(pygmds.DIM3 | pygmds.F | pygmds.E |
                       pygmds.N | pygmds.F2N | pygmds.F2E |
                       pygmds.F2N | pygmds.E2N | pygmds.N2F)
m = pygmds.Mesh(mod)
ios = pygmds.IOService(m)

# File reading
reader = pygmds.VTKReader(ios)
reader.setCellOptions(pygmds.MeshModel(pygmds.N | pygmds.F))
reader.read(file_name)
if m.getNbNodes() != 222:
    sys.exit(1)
if m.getNbFaces() != 398:
    sys.exit(2)

#in success, return 0 for ctest
sys.exit(0)
