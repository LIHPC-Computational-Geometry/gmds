import sys
from gmds import *


def setUp(fac_mng, file_name):
    model_3d = MeshModel(DIM3 | R| F | E | N | R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E)
    m_vol = Mesh(model_3d)

    ios = IOService(m_vol)
    # File reading
    reader = VTKReader(ios)
    reader.set_cell_options(MeshModel(N|R))
    reader.read(file_name)

    doc = MeshDoctor(m_vol)
    doc.build_faces_and_R2F()
    doc.build_edges_and_X2E()
    doc.update_upward_connectivity()

    fac_mng.init_from_3d_mesh(m_vol)

def test_blocking():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    blocking = Blocking(fm,True)
    blocks = blocking.get_all_blocks()
    # We've got a single block which is the bounding box of the model
    assert len(blocks) == 1
    assert blocking.get_nb_blocks() == 1
    assert blocking.get_nb_faces() == 6
    assert blocking.get_nb_edges() == 12
    assert blocking.get_nb_nodes() == 8

    edges = blocking.get_all_edges()
    blocking.cut_sheet(edges[0])
    assert blocking.get_nb_blocks() == 2
    assert blocking.get_nb_faces() == 11
    assert blocking.get_nb_edges() == 20
    assert blocking.get_nb_nodes() == 12

    blocking.remove_block(blocking.get_all_blocks()[1])
    assert blocking.get_nb_blocks() == 1
    assert blocking.get_nb_faces() == 6
    assert blocking.get_nb_edges() == 12
    assert blocking.get_nb_nodes() == 8

    mesh_model = MeshModel(DIM3 | R| F | E | N | R2N |  F2N | E2N)
    mesh_block = Mesh(mesh_model)
    blocking.convert_to_mesh(mesh_block)
    assert  mesh_block.get_nb_regions() == 1
    assert  mesh_block.get_nb_faces()   == 6
    assert  mesh_block.get_nb_edges()   == 12
    assert  mesh_block.get_nb_nodes()   == 8


def test_blocking_with_classification():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    blocking = Blocking(fm,True)
    classifier = BlockingClassifier(blocking)
    errors = classifier.classify(0.01,0.1)

    assert len(errors.non_captured_points) == 0
    assert len(errors.non_captured_curves) == 0
    assert len(errors.non_captured_surfaces) == 0
