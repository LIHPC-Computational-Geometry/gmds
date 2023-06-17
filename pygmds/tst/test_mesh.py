import sys
import gmds

def test_import_vtk():
    # the second argument is the directory where test samples are put
    file_name = str(sys.argv[2]) + "/quarter_crown.vtk"
    print("file  :" + file_name)
    mod = gmds.MeshModel(gmds.DIM3 | gmds.F | gmds.E |
                           gmds.N | gmds.F2N | gmds.F2E |
                           gmds.F2N | gmds.E2N | gmds.N2F)
    m = gmds.Mesh(mod)
    ios = gmds.IOService(m)

    # File reading
    reader = gmds.VTKReader(ios)
    reader.set_cell_options(gmds.MeshModel(gmds.N | gmds.F))
    reader.read(file_name)
    assert m.get_nb_nodes() == 222
    assert m.get_nb_faces() == 398
