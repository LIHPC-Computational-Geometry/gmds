import sys
import gmds.mesh

def test_import_vtk():
    # the second argument is the directory where test samples are put
    file_name = str(sys.argv[2]) + "/quarter_crown.vtk"
    print("file  :" + file_name)
    mod = gmds.mesh.MeshModel(gmds.mesh.DIM3 | gmds.mesh.F | gmds.mesh.E |
                           gmds.mesh.N | gmds.mesh.F2N | gmds.mesh.F2E |
                           gmds.mesh.F2N | gmds.mesh.E2N | gmds.mesh.N2F)
    m = gmds.mesh.Mesh(mod)
    ios = gmds.mesh.IOService(m)

    # File reading
    reader = gmds.mesh.VTKReader(ios)
    reader.set_cell_options(gmds.mesh.MeshModel(gmds.mesh.N | gmds.mesh.F))
    reader.read(file_name)
    assert m.get_nb_nodes() == 222
    assert m.get_nb_faces() == 398
