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

def test_convert_blocking_into_mesh():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    blocking = Blocking(fm,True)
    classifier = BlockingClassifier(blocking)
    errors = classifier.classify(0.01,0.1)
    m=Mesh(MeshModel(DIM3|R|F|E|N|F2N|R2N|E2N))
    blocking.convert_to_mesh(m)
    assert m.get_nb_nodes() == 8
    assert m.get_nb_faces() == 6
def test_sheet_edges():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    blocking = Blocking(fm,True)
    sheet_sets = blocking.get_all_sheet_edge_sets()

    assert len(sheet_sets) == 3
    for ss in sheet_sets:
        assert len(ss) == 4

    p = Point(5,5,5)
    for  sh_edges in sheet_sets:
        assert len(sh_edges) == 4
        dist_coord = blocking.get_projection_info(p,sh_edges)
        min_dist = 1000.0
        for dc in dist_coord:
            if dc[0]<min_dist:
                 min_dist = dc[0]

        assert min_dist <0.01


def test_all_ids():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    bl = Blocking(fm,True)

    classifier = BlockingClassifier(bl)
    errors = classifier.classify(0.01,0.1)

    node_ids = bl.get_all_id_nodes()
    for id in node_ids:
        geom_dim, geom_id, p = bl.get_node_info(id)
        # all the node are classified on geom point
        assert geom_dim == 0
        assert abs(p.x()) == 5
        assert abs(p.y()) == 5
        assert abs(p.z()) == 5


def test_face_center_normal():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    blocking = Blocking(fm,True)
    assert blocking.get_nb_faces() == 6
    faces = blocking.get_all_faces()
    assert len(faces) == 6
    for f in faces:
        n = blocking.get_normal_of_face(f)
        assert (abs(n.x()) < 0.001) or (abs(n.x())-1 < 0.001)
        assert (abs(n.y()) < 0.001) or (abs(n.y())-1 < 0.001)
        assert (abs(n.z()) < 0.001) or (abs(n.z())-1 < 0.001)
