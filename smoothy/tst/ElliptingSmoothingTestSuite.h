/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/smoothy/EllipticSmoothing.h>
#include <gmds/smoothy/EllipticSmoother2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
#include <random>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
/*----------------------------------------------------------------------------*/
TEST(EllipticSmoothingTestSuite, grid2D_internal_smoothing)
{
	Mesh m(MeshModel(DIM3 | F | E | N | F2N | F2E | E2F | E2N | N2E | N2F));

	GridBuilder gb(&m,2);
	gb.execute(5,1,5,1);

	MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	doc.orient2DFaces();

	for(auto f_id:m.faces()) {
		Face f=m.get<Face>(f_id);
		if (f.normal().dot(math::Vector3d({.0, .0, 1.0})) <= 0) {
			std::vector<TCellID> ns = f.getIDs<Node>();
			std::vector<TCellID> ns2(ns.size());
			for (auto i = 0; i < ns.size(); i++)
				ns2[ns.size() - 1 - i] = ns[i];
			f.set<Node>(ns2);
		}
	}
	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator2D op(&m);
	auto mark_node_NAN = m.newMark<Node>();
	auto mark_node_on_pnt = m.newMark<Node>();
	auto mark_node_on_crv = m.newMark<Node>();
	auto mark_edge_on_crv = m.newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	auto nb_locked = 0;
	auto mark_bnd_nodes = m.newMark<Node>();
	for (auto n_id : m.nodes()) {
		if (m.isMarked<Node>(n_id, mark_node_on_crv) || m.isMarked<Node>(n_id, mark_node_on_pnt)) {
			nb_locked += 1;
			m.mark<Node>(n_id, mark_bnd_nodes);
		}
	}
	//==================================================================
	// PERFORM THE PERTURBATION
	//==================================================================
	constexpr int FLOAT_MIN = -100;
	constexpr int FLOAT_MAX = 108;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);

	for(auto n_id:m.nodes()){
		if(!m.isMarked<Node>(n_id,mark_bnd_nodes)) {
			Node n = m.get<Node>(n_id);
			n.setXYZ(distr(eng), distr(eng), 0);
		}
	}
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================

	std::vector<double> node_coords(2*m.getNbNodes());
	std::vector<bool> lock_nodes(2*m.getNbNodes());
	for(auto n_id:m.nodes()){
		math::Point pi = m.get<Node>(n_id).point();
		node_coords[2*n_id]=pi.X();
		node_coords[2*n_id+1]=pi.Y();
		lock_nodes[2*n_id]= (m.isMarked<Node>(n_id, mark_bnd_nodes));
		lock_nodes[2*n_id+1]= (m.isMarked<Node>(n_id, mark_bnd_nodes));
	}

	std::vector<std::array<int, 3>> triangles(4*m.getNbFaces());
	constexpr int quadrature[4][3]= {{0,1,3},{1,2,0},{2,3,1},{3,0,2}};
	for(auto f_id:m.faces()){
		std::vector<TCellID> nids = m.get<Face>(f_id).getIDs<Node>();
		for(auto i=0;i<4;i++) {
			triangles[4 * f_id + i] = {static_cast<int>(nids[quadrature[i][0]]),
			                           static_cast<int>(nids[quadrature[i][1]]),
			                           static_cast<int>(nids[quadrature[i][2]])};
		}
	}
	EllipticSmoothingMeshGlue2D tri(node_coords,triangles,lock_nodes);

	EllipticSmoothingOptions options = ESO_default2D;
	options.maxiter=10;//00;
	options.eps_from_theorem=true;
	options.theta = 1e-3;
	options.bfgs_threshold = 1e-9;
	options.debug=0;
	EllipticSmoothing2D smoother2D(tri, triangles.size(),options);
	//	smoother2D.start_eps = 1E-4;
	smoother2D.execute();
	tri.get_verts(node_coords);
	for(auto n_id:m.nodes()){
		m.get<Node>(n_id).setXYZ(node_coords[2*n_id],node_coords[2*n_id+1],.0);
	}
	for(auto i:m.nodes()){
		Node ni = m.get<Node>(i);
		ASSERT_TRUE(ni.X()>=0. && ni.X()<=4.);
		ASSERT_TRUE(ni.Y()>=0. && ni.Y()<=4.);
		ASSERT_TRUE(ni.Z()>=0. && ni.Z()<=4.);
	}
	for(auto j:m.faces()){
		Face fj = m.get<Face>(j);
		ASSERT_NEAR(1.0,fj.computeScaledJacobian2D(),0.01);
	}
}
/*----------------------------------------------------------------------------*/
TEST(EllipticSmoothingTestSuite, grid2D_smoother)
{
	Mesh m(MeshModel(DIM3 | F | E | N | F2N | F2E | E2F | E2N | N2E | N2F));

	GridBuilder gb(&m,2);
	gb.execute(5,1,5,1);

	MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();
	for(auto f_id:m.faces()) {
		Face f=m.get<Face>(f_id);
		if (f.normal().dot(math::Vector3d({.0, .0, 1.0})) <= 0) {
			std::vector<TCellID> ns = f.getIDs<Node>();
			std::vector<TCellID> ns2(ns.size());
			for (auto i = 0; i < ns.size(); i++)
				ns2[ns.size() - 1 - i] = ns[i];
			f.set<Node>(ns2);
		}
	}

	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator2D op(&m);
	auto mark_node_NAN = m.newMark<Node>();
	auto mark_node_on_pnt = m.newMark<Node>();
	auto mark_node_on_crv = m.newMark<Node>();
	auto mark_edge_on_crv = m.newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	auto nb_locked = 0;
	auto mark_bnd_nodes = m.newMark<Node>();
	for (auto n_id : m.nodes()) {
		if (m.isMarked<Node>(n_id, mark_node_on_crv) || m.isMarked<Node>(n_id, mark_node_on_pnt)) {
			nb_locked += 1;
			m.mark<Node>(n_id, mark_bnd_nodes);
		}
	}
	//==================================================================
	// PERFORM THE PERTURBATION
	//==================================================================
	constexpr int FLOAT_MIN = -100;
	constexpr int FLOAT_MAX = 108;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);

	for(auto n_id:m.nodes()){
		if(!m.isMarked<Node>(n_id,mark_bnd_nodes)) {
			Node n = m.get<Node>(n_id);
			n.setXYZ(distr(eng), distr(eng), 0);
		}
	}
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOWEllipticSmoothingTestSuite.grid2D_smoother (1)
	//==================================================================

	EllipticSmoother2D smoother2D(&m);
	smoother2D.lock(mark_bnd_nodes);
	smoother2D.execute();

	for(auto i:m.nodes()){
		Node ni = m.get<Node>(i);
		ASSERT_TRUE(ni.X()>=0. && ni.X()<=4.);
		ASSERT_TRUE(ni.Y()>=0. && ni.Y()<=4.);
		ASSERT_TRUE(ni.Z()>=0. && ni.Z()<=4.);
	}
	for(auto j:m.faces()){
		Face fj = m.get<Face>(j);
		ASSERT_NEAR(1.0,fj.computeScaledJacobian2D(),0.01);
	}
}
/*----------------------------------------------------------------------------*/
TEST(EllipticSmoothingTestSuite, grid2D_smoother_UC)
{
	Mesh m(MeshModel(DIM3 | F | E | N | F2N | F2E | E2F | E2N | N2E | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1,3,1);

	MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();
	for(auto f_id:m.faces()) {
		Face f=m.get<Face>(f_id);
		if (f.normal().dot(math::Vector3d({.0, .0, 1.0})) <= 0) {
			std::vector<TCellID> ns = f.getIDs<Node>();
			std::vector<TCellID> ns2(ns.size());
			for (auto i = 0; i < ns.size(); i++)
				ns2[ns.size() - 1 - i] = ns[i];
			f.set<Node>(ns2);
		}
	}

	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator2D op(&m);
	auto mark_node_NAN = m.newMark<Node>();
	auto mark_node_on_pnt = m.newMark<Node>();
	auto mark_node_on_crv = m.newMark<Node>();
	auto mark_edge_on_crv = m.newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	auto nb_locked = 0;
	auto mark_bnd_nodes = m.newMark<Node>();
	for (auto n_id : m.nodes()) {
		if (m.isMarked<Node>(n_id, mark_node_on_crv) || m.isMarked<Node>(n_id, mark_node_on_pnt)) {
			nb_locked += 1;
			m.mark<Node>(n_id, mark_bnd_nodes);
		}
	}
	//==================================================================
	// PERFORM THE PERTURBATION
	//==================================================================
	constexpr int FLOAT_MIN = -100;
	constexpr int FLOAT_MAX = 108;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);

	for(auto n_id:m.nodes()){
		if(!m.isMarked<Node>(n_id,mark_bnd_nodes)) {
			Node n = m.get<Node>(n_id);
			n.setXYZ(distr(eng), distr(eng), 0);
		}
	}
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOWEllipticSmoothingTestSuite.grid2D_smoother (1)
	//==================================================================
	EllipticSmoother2D smoother2D(&m);
	smoother2D.lock(mark_bnd_nodes);
	smoother2D.execute();

	for(auto i:m.nodes()){
		Node ni = m.get<Node>(i);
		ASSERT_TRUE(ni.X()>=0. && ni.X()<=4.);
		ASSERT_TRUE(ni.Y()>=0. && ni.Y()<=4.);
		ASSERT_TRUE(ni.Z()>=0. && ni.Z()<=4.);
	}
	for(auto j:m.faces()){
		Face fj = m.get<Face>(j);
		ASSERT_NEAR(1.0,fj.computeScaledJacobian2D(),0.01);
	}


}

/*----------------------------------------------------------------------------*/
TEST(EllipticSmoothingTestSuite, aero_test1)
{
	Mesh m(MeshModel(DIM3 | F | E | N | F2N | F2E | E2F | E2N | N2E | N2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Aero/2D/ElipticSmoothTest1.vtk";
	IGMeshIOService ioService(&m);

	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();
	for(auto f_id:m.faces()) {
		Face f=m.get<Face>(f_id);
		if (f.normal().dot(math::Vector3d({.0, .0, 1.0})) <= 0) {
			std::vector<TCellID> ns = f.getIDs<Node>();
			std::vector<TCellID> ns2(ns.size());
			for (auto i = 0; i < ns.size(); i++)
				ns2[ns.size() - 1 - i] = ns[i];
			f.set<Node>(ns2);
		}
	}

	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator2D op(&m);
	auto mark_node_NAN = m.newMark<Node>();
	auto mark_node_on_pnt = m.newMark<Node>();
	auto mark_node_on_crv = m.newMark<Node>();
	auto mark_edge_on_crv = m.newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	auto nb_locked = 0;
	auto mark_bnd_nodes = m.newMark<Node>();
	for (auto n_id : m.nodes()) {
		if (m.isMarked<Node>(n_id, mark_node_on_crv) || m.isMarked<Node>(n_id, mark_node_on_pnt)) {
			nb_locked += 1;
			m.mark<Node>(n_id, mark_bnd_nodes);
		}
	}
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================
	EllipticSmoother2D smoother2D(&m);
	smoother2D.lock(mark_bnd_nodes);
	smoother2D.execute();

	/*
	   VTKWriter w(&ioService);
	   w.setCellOptions(gmds::N|gmds::F);
	   w.setDataOptions(gmds::N|gmds::F);
	   w.write("out_aero2.vtk");
	*/
}

/*----------------------------------------------------------------------------*/
TEST(EllipticSmoothingTestSuite, aero_test2)
{
	Mesh m(MeshModel(DIM3 | F | E | N | F2N | F2E | E2F | E2N | N2E | N2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Aero/2D/ElipticSmoothTest2.vtk";
	IGMeshIOService ioService(&m);

	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();
	for(auto f_id:m.faces()) {
		Face f=m.get<Face>(f_id);
		if (f.normal().dot(math::Vector3d({.0, .0, 1.0})) <= 0) {
			std::vector<TCellID> ns = f.getIDs<Node>();
			std::vector<TCellID> ns2(ns.size());
			for (auto i = 0; i < ns.size(); i++)
				ns2[ns.size() - 1 - i] = ns[i];
			f.set<Node>(ns2);
		}
	}

	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator2D op(&m);
	auto mark_node_NAN = m.newMark<Node>();
	auto mark_node_on_pnt = m.newMark<Node>();
	auto mark_node_on_crv = m.newMark<Node>();
	auto mark_edge_on_crv = m.newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	auto nb_locked = 0;
	auto mark_bnd_nodes = m.newMark<Node>();
	for (auto n_id : m.nodes()) {
		if (m.isMarked<Node>(n_id, mark_node_on_crv) || m.isMarked<Node>(n_id, mark_node_on_pnt)) {
			nb_locked += 1;
			m.mark<Node>(n_id, mark_bnd_nodes);
		}
	}
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================
	EllipticSmoother2D smoother2D(&m);
	smoother2D.lock(mark_bnd_nodes);
	smoother2D.execute();

	/*
	   VTKWriter w(&ioService);
	   w.setCellOptions(gmds::N|gmds::F);
	   w.setDataOptions(gmds::N|gmds::F);
	   w.write("out_aero2.vtk");
	*/
}

