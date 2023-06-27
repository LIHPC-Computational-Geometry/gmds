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

def test_point_queries():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    assert fm.get_nb_points() == 8
    pnts = fm.get_points()
    for p in pnts:
        assert abs(p.x()) == 5
        assert abs(p.y()) == 5
        assert abs(p.z()) == 5
        assert len(p.curves()) == 3
        assert len(p.surfaces()) == 3
        assert len(p.volumes()) == 1
def test_curve_queries():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    assert fm.get_nb_curves() == 12
    crvs =  fm.get_curves()
    assert len(crvs) == 12
    for c in crvs:
        assert c.length() ==10
        assert c.curvature_info() == CURVATURE_CONVEX
        assert len(c.points()) == 2
        assert len(c.surfaces()) == 2
        assert len(c.volumes()) == 1
        tangent_0 = c.tangent(0)
        tangent_1 = c.tangent(1)
        p0 = c.points()[0]
        p1 = c.points()[1]
        v01 = Vector3d(p1.x()-p0.x(),p1.y()-p0.y(),p1.z()-p0.z())
        v01 = v01.get_normalize()
        assert abs(tangent_0.dot(v01))-1<0.001
        assert abs(tangent_1.dot(v01))-1<0.001
def test_surface_queries():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    assert fm.get_nb_surfaces() == 6
    surfs =  fm.get_surfaces()
    assert len(surfs) == 6
    for s in surfs:
        assert abs(s.area() - 100) < 0.001
        assert len(s.points()) == 4
        assert len(s.curves()) == 4
        assert len(s.volumes()) == 1
def test_volume_queries():
    input_geom_file = str(sys.argv[2]) + "/tet_in_box.vtk"
    fm = FACManager()
    setUp(fm, input_geom_file)
    assert fm.get_nb_volumes() == 1
    vols =  fm.get_volumes()
    assert len(vols) == 1
    for v in vols:
        assert len(v.points()) == 8
        assert len(v.curves()) == 12
        assert len(v.surfaces()) == 6
        bb  = v.bbox()
        assert bb[0] == -5
        assert bb[1] == -5
        assert bb[2] == -5
        assert bb[3] ==  5
        assert bb[4] ==  5
        assert bb[5] ==  5