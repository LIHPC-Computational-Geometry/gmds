//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/smoothy/LaplacianSmoother.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(GeomSmootherTestSuite, tet_in_cube)
{
	// WE WRITE
	Mesh m_vol(gmds::MeshModel(DIM3 | R | F | E | N | R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | N2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/tet_in_box.vtk";
	IGMeshIOService ioService(&m_vol);

	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);

	MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	cad::FACManager manager;
	cad::GeomMeshLinker linker;

	manager.initAndLinkFrom3DMesh(&m_vol, &linker);

	smoothy::LaplacianSmoother smoother(&m_vol, &linker);
	smoother.setNbIterations(100);
	smoother.setNodesAll();

	// We perturb mesh node locations to test the smoothing algorithm
	for (auto n_id : m_vol.nodes()) {
		if (n_id == 8) {
			m_vol.get<Node>(n_id).setPoint(math::Point(5, 1, 5));
		}
		else if (n_id == 10) {
			m_vol.get<Node>(n_id).setPoint(math::Point(-5, 1, 5));
		}
		else if (n_id == 20) {
			m_vol.get<Node>(n_id).setPoint(math::Point(-2, 2, 5));
		}
		else if (n_id == 28) {
			m_vol.get<Node>(n_id).setPoint(math::Point(1, 1, 3));
		}
	}

	std::string vtk_file2 = ("TOTO_moved.vtk");

	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::R);
	vtkWriter.setDataOptions(gmds::N | gmds::R);
	vtkWriter.write(vtk_file2);

	smoother.smoothCurves();

	std::string vtk_file3 = ("TOTO_smooth.vtk");

	vtkWriter.write(vtk_file3);

	ASSERT_NEAR(m_vol.get<Node>(8).point().Y(), 0, 0.01);
	ASSERT_NEAR(m_vol.get<Node>(10).point().Y(), 0, 0.01);

	smoother.setNbIterations(10);
	smoother.smoothSurfaces();

	ASSERT_NEAR(m_vol.get<Node>(20).point().X(), -1, 0.01);
	ASSERT_NEAR(m_vol.get<Node>(20).point().Y(), 1, 0.01);

	smoother.setNbIterations(1);
	smoother.smoothVolumes();

	ASSERT_NEAR(m_vol.get<Node>(28).point().X(), 2.08, 0.01);
	ASSERT_NEAR(m_vol.get<Node>(28).point().Y(), 0.64, 0.01);
	ASSERT_NEAR(m_vol.get<Node>(28).point().Z(), 1.97, 0.01);
}

/*----------------------------------------------------------------------------*/
TEST(GeomSmootherTestSuite, test3)
{
	Mesh m_vol(gmds::MeshModel(DIM3 | R | F | E | N | R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | N2F));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/hexa.vtk";

	IGMeshIOService ioService(&m_vol);

	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);

	m_vol.deleteRegion(0);
	m_vol.deleteRegion(1);
	m_vol.deleteRegion(2);
	m_vol.deleteRegion(3);

	MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	cad::FACManager manager;
	cad::GeomMeshLinker linker;

	manager.initAndLinkFrom3DMesh(&m_vol, &linker);

	smoothy::LaplacianSmoother smoother(&m_vol,&linker);
smoother.setNbIterations(10);
smoother.setNodesAll();
	smoother.smoothCurves();

	smoother.smoothSurfaces();

	smoother.smoothVolumes();

	std::string vtk_file2 = ("test_samples/hexa_out.vtk");

	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::R);
	// vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write(vtk_file2);
	ASSERT_TRUE(true);
}

/*----------------------------------------------------------------------------*/
TEST(GeomSmootherTestSuite, DISABLED_test_notch_refined)
{
	std::string dir(TEST_SAMPLES_DIR);
	std::string file_geom = dir + "/Notch/notch_tet.vtk";
	std::string file_mesh = dir + "/Notch/notch_hexa.vtk";

	auto nb_iterations = 10;
	//==================================================================
	// GEOMETRY READING
	//==================================================================
	Mesh geometry(MeshModel(DIM3 | R | F | E | N | R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | N2F));

	IGMeshIOService ioService(&geometry);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N | R);
	vtkReader.read(file_geom);
	MeshDoctor doc(&geometry);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	cad::FACManager geom_manager;
	geom_manager.initFrom3DMesh(&geometry);
	//==================================================================
	// MESH READING
	//==================================================================
	// the used model is specified according to the geom smoother requirements.
	Mesh m(MeshModel(DIM3 | R | F | E | N | R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | N2F));

	IGMeshIOService ioService2(&m);
	VTKReader vtkReader2(&ioService2);
	vtkReader2.setCellOptions(N | R);
	vtkReader2.read(file_mesh);
	MeshDoctor doc2(&m);
	doc2.buildFacesAndR2F();
	doc2.buildEdgesAndX2E();
	doc2.updateUpwardConnectivity();

	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator op(&m);
	auto mark_node_NAN = m.newMark<Node>();
	auto mark_node_on_pnt = m.newMark<Node>();
	auto mark_node_on_crv = m.newMark<Node>();
	auto mark_node_on_srf = m.newMark<Node>();
	auto mark_edge_on_crv = m.newMark<Edge>();
	auto mark_edge_on_srf = m.newMark<Edge>();
	auto mark_face_on_srf = m.newMark<Face>();

	op.markCellOnGeometry(mark_face_on_srf, mark_edge_on_srf, mark_node_on_srf, mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);

	//==================================================================
	// TRY TO CLASSIFY MESH CELLS
	// The classification process is not straightforward as we have to
	// consider that the mesh to classify is a "coarse" version of the
	// geometry. For instance, if you try and project the middle point
	// of an edge onto a curve, the distance between the initial point
	// and the projected one will be probably not equal to zero. While,
	// the same point will be lying on adjacent surface with a null
	// distance. So the strategy is to
	// 1. Classify faces on surfaces
	// 2. Get edges that are shared by two faces classified onto 2
	//    surfaces. Those edges must be classified on curves. It is
	//    also the case of the edge endpoint
	// 3. Do the same for nodes shared by at least 3 faces classified
	//    onto 3 surfaces. Those nodes must be classified onto points.
	//==================================================================

	cad::GeomMeshLinker linker(&m, &geom_manager);
	std::vector<cad::GeomSurface *> surfaces;
	std::vector<cad::GeomCurve *> curves;
	std::vector<cad::GeomPoint *> points;

	geom_manager.getSurfaces(surfaces);
	geom_manager.getCurves(curves);
	geom_manager.getPoints(points);

	//==================================================================
	// First, we classify each face
	for (auto f_id : m.faces()) {
		Face f = m.get<Face>(f_id);
		if (m.isMarked(f, mark_face_on_srf)) {
			// we've got a boundary face
			math::Point p = f.center();

			double min_dist = 100000;
			int min_entity_dim = -1;
			int min_entity_id = -1;
			for (auto s : surfaces) {
				math::Point closest_pnt = p;
				s->project(closest_pnt);
				double dist = p.distance2(closest_pnt);
				if (dist < min_dist) {
					min_dist = dist;
					min_entity_dim = 2;
					min_entity_id = s->id();
				}
			}
			if (min_entity_dim == 2) {
				linker.linkFaceToSurface(f_id, min_entity_id);
			}
		}
	}
	//==================================================================
	// Second, we classify each edge
	for (auto e_id : m.edges()) {
		Edge e = m.get<Edge>(e_id);
		if (m.isMarked(e, mark_edge_on_crv) || m.isMarked(e, mark_edge_on_srf)) {
			// we've got a boundary edge, now we get the 2 boundary faces
			//  around e
			std::vector<Node> e_nodes = e.get<Node>();
			std::vector<TCellID> adj_face_id = e.getIDs<Face>();
			std::vector<TCellID> adj_bnd_face_ids;
			for (auto f_id : adj_face_id) {
				if (m.isMarked<Face>(f_id, mark_face_on_srf)) {
					adj_bnd_face_ids.push_back(f_id);
				}
			}
			auto surf0_id = linker.getGeomId<Face>(adj_bnd_face_ids[0]);
			auto surf1_id = linker.getGeomId<Face>(adj_bnd_face_ids[1]);

			if (surf0_id == surf1_id) {
				// edge embedded in the surface
				linker.linkEdgeToSurface(e_id, surf1_id);
			}
			else {
				// it must be an edge classified on curves
				// we look for the closest curve now
				math::Point p0 = e_nodes[0].point();
				math::Point p1 = e_nodes[1].point();
				math::Point p = 0.5 * (p0 + p1);

				double min_dist = 100000;
				int min_entity_dim = -1;
				int min_entity_id = -1;
				for (auto c : curves) {
					math::Point closest_pnt = p;
					c->project(closest_pnt);
					double dist = p.distance2(closest_pnt);
					if (dist < min_dist) {
						min_dist = dist;
						min_entity_dim = 1;
						min_entity_id = c->id();
					}
				}
				linker.linkEdgeToCurve(e_id, min_entity_id);
			}
		}
	}
	//==================================================================
	// we classify each node
	auto on_pnt = 0, on_curve = 0, on_surf = 0;
	for (auto n_id : m.nodes()) {
		Node n = m.get<Node>(n_id);
		if (m.isMarked(n, mark_node_on_pnt) || m.isMarked(n, mark_node_on_crv) || m.isMarked(n, mark_node_on_srf)) {
			// we've got a boundary node
			math::Point node_loc = n.point();
			double min_dist = 100000;
			int min_entity_dim = -1;
			int min_entity_id = -1;
			for (auto p : points) {
				math::Point closest_pnt = node_loc;
				double dist = node_loc.distance(p->point());
				if (dist < min_dist) {
					min_dist = dist;
					min_entity_dim = 0;
					min_entity_id = p->id();
				}
			}
			for (auto c : curves) {
				math::Point closest_pnt = node_loc;
				c->project(closest_pnt);
				double dist = node_loc.distance(closest_pnt);
				if (dist < min_dist) {
					// WARNING: Take care of this trick that is not good at all but mandatory to be staying on the
					//  curve and not on the surface
					if (dist < 1e-4) dist = 0;
					min_dist = dist;
					min_entity_dim = 1;
					min_entity_id = c->id();
				}
			}
			for (auto s : surfaces) {
				math::Point closest_pnt = node_loc;
				s->project(closest_pnt);
				double dist = node_loc.distance(closest_pnt);
				if (dist < min_dist) {
					min_dist = dist;
					min_entity_dim = 2;
					min_entity_id = s->id();
				}
			}

			if (min_entity_dim == 0) {
				on_pnt++;
				linker.linkNodeToPoint(n_id, min_entity_id);
			}
			else if (min_entity_dim == 1) {
				on_curve++;
				linker.linkNodeToCurve(n_id, min_entity_id);
			}
			else if (min_entity_dim == 2) {
				on_surf++;
				linker.linkNodeToSurface(n_id, min_entity_id);
			}
			else {
				throw GMDSException("Link error for classifying a node");
			}
		}
	}
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================
	smoothy::LaplacianSmoother smoother(&m,&linker);
	smoother.setNbIterations(nb_iterations);
	smoother.setNodesAll();
	smoother.smoothCurves();
	smoother.smoothSurfaces();
	smoother.smoothVolumes();

	ASSERT_TRUE(true);
}
