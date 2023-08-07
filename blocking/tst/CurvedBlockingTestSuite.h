/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
void
setUp(gmds::cad::FACManager &AGeomManager)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N |
	                                 gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N |
	                                 gmds::F2R | gmds::F2E
	                                 | gmds::E2F | gmds::E2N | gmds::N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/tet_in_box.vtk";
	gmds::IGMeshIOService ioService(&m_vol);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	gmds::MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	AGeomManager.initFrom3DMesh(&m_vol);
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, global_cell_accessors)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model);
	gmds::math::Point p000(0, 0, 0);
	gmds::math::Point p010(0, 1, 0);
	gmds::math::Point p110(1, 1, 0);
	gmds::math::Point p100(1, 0, 0);

	gmds::math::Point p001(0, 0, 1);
	gmds::math::Point p011(0, 1, 1);
	gmds::math::Point p111(1, 1, 1);
	gmds::math::Point p101(1, 0, 1);

	gmds::math::Point p002(0, 0, 2);
	gmds::math::Point p012(0, 1, 2);
	gmds::math::Point p112(1, 1, 2);
	gmds::math::Point p102(1, 0, 2);
	auto b1 = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);
	auto b2 = bl.create_block(p001, p011, p111, p101, p002, p012, p112, p102);
	ASSERT_EQ(16, bl.get_nb_cells<0>());
	ASSERT_EQ(24, bl.get_nb_cells<1>());
	ASSERT_EQ(12, bl.get_nb_cells<2>());
	ASSERT_EQ(2,  bl.get_nb_cells<3>());
	ASSERT_EQ(16, bl.get_all_nodes().size());
	ASSERT_EQ(24, bl.get_all_edges().size());
	ASSERT_EQ(12, bl.get_all_faces().size());
	ASSERT_EQ(2,  bl.get_all_blocks().size());
	bl.sew<3>(b1->dart(), b2->dart());
	ASSERT_EQ(12, bl.get_nb_cells<0>());
	ASSERT_EQ(20, bl.get_nb_cells<1>());
	ASSERT_EQ(11, bl.get_nb_cells<2>());
	ASSERT_EQ(2, bl.get_nb_cells<3>());
	ASSERT_EQ(12, bl.get_all_nodes().size());
	ASSERT_EQ(20, bl.get_all_edges().size());
	ASSERT_EQ(11, bl.get_all_faces().size());
	ASSERT_EQ(2,  bl.get_all_blocks().size());
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, remove_block)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model);
	gmds::math::Point p000(0, 0, 0);
	gmds::math::Point p010(0, 1, 0);
	gmds::math::Point p110(1, 1, 0);
	gmds::math::Point p100(1, 0, 0);

	gmds::math::Point p001(0, 0, 1);
	gmds::math::Point p011(0, 1, 1);
	gmds::math::Point p111(1, 1, 1);
	gmds::math::Point p101(1, 0, 1);

	gmds::math::Point p002(0, 0, 2);
	gmds::math::Point p012(0, 1, 2);
	gmds::math::Point p112(1, 1, 2);
	gmds::math::Point p102(1, 0, 2);
	auto b1 = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);
	auto b2 = bl.create_block(p001, p011, p111, p101, p002, p012, p112, p102);
	bl.sew<3>(b1->dart(), b2->dart());
	ASSERT_EQ(bl.get_nb_cells<3>(),2);
	bl.remove_block(b1);
	ASSERT_EQ(bl.get_nb_cells<3>(),1);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, single_block)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model);
	gmds::math::Point p000(0, 0, 0);
	gmds::math::Point p010(0, 1, 0);
	gmds::math::Point p110(1, 1, 0);
	gmds::math::Point p100(1, 0, 0);

	gmds::math::Point p001(0, 0, 1);
	gmds::math::Point p011(0, 1, 1);
	gmds::math::Point p111(1, 1, 1);
	gmds::math::Point p101(1, 0, 1);

	gmds::math::Point p002(0, 0, 2);
	gmds::math::Point p012(0, 1, 2);
	gmds::math::Point p112(1, 1, 2);
	gmds::math::Point p102(1, 0, 2);
	auto b = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);
	ASSERT_EQ(b->info().geom_dim, 4);
	ASSERT_EQ(b->info().geom_id, -1);
	auto fs = bl.get_faces_of_block(b);
	ASSERT_EQ(fs.size(), 6);
	auto block_center = bl.get_center_of_block(b);
	for (auto f : fs) {
		gmds::math::Point face_center = bl.get_center_of_face(f);
		ASSERT_NEAR(block_center.distance(face_center), 0.5, 1e-8);
	}
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, single_block_parallel_edges)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model);
	gmds::math::Point p000(0, 0, 0);
	gmds::math::Point p010(0, 1, 0);
	gmds::math::Point p110(1, 1, 0);
	gmds::math::Point p100(1, 0, 0);

	gmds::math::Point p001(0, 0, 1);
	gmds::math::Point p011(0, 1, 1);
	gmds::math::Point p111(1, 1, 1);
	gmds::math::Point p101(1, 0, 1);

	bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);

	for (auto it = bl.gmap()->attributes<1>().begin(), itend = bl.gmap()->attributes<1>().end(); it != itend; ++it) {
		std::vector<gmds::blocking::CurvedBlocking::Edge> parallel_edges;
		bl.get_all_sheet_edges(it,parallel_edges);
		ASSERT_EQ(4, parallel_edges.size());
	}

	std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge> > all_edges= bl.get_all_sheet_edge_sets();
	ASSERT_EQ(3, all_edges.size());
	for(auto sh_edges: all_edges){
		ASSERT_EQ(4, sh_edges.size());
	}
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, get_edges_of_a_block)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model, true);
	auto b = bl.get_all_blocks()[0];
	auto edges = bl.get_edges_of_block(b);
	ASSERT_EQ(12,edges.size());
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, split_one_block_twice)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model, true);
	gmds::blocking::CurvedBlockingClassifier cl(&bl);
	cl.classify();

	auto e = bl.get_all_edges()[0];
	auto e_id = e->info().geom_id;
	auto e_dim = e->info().geom_dim;

	auto e2 = bl.gmap()->attribute<1>(bl.gmap()->alpha<1>(e->dart()));
	bl.cut_sheet(e);
	ASSERT_EQ(12,bl.get_nb_cells<0>());
	ASSERT_EQ(20,bl.get_nb_cells<1>());
	ASSERT_EQ(11,bl.get_nb_cells<2>());
	ASSERT_EQ(2,bl.get_nb_cells<3>());
	ASSERT_TRUE(bl.gmap()->is_valid());

	//after splitting e, one node and 2 edges must be classified on its original geometric cell
	auto classified_nodes = 0;
	auto classified_edges = 0;
	for(auto cur_edge:bl.get_all_edges()){
		if(cur_edge->info().geom_id==e_id && cur_edge->info().geom_dim==e_dim)
			classified_edges++;
	}
	for(auto cur_node:bl.get_all_nodes()){
		if(cur_node->info().geom_id==e_id && cur_node->info().geom_dim==e_dim)
			classified_nodes++;
	}
	ASSERT_EQ(2, classified_edges);
	ASSERT_EQ(1, classified_nodes);
	//we check the attribute values
	gmds::math::Point p_cut(-5,5,2);
	bl.cut_sheet(e2, p_cut);
	//We cut along Z axis and point are so located on Z-plane with Z=5, -5, or 2.
	for(auto n: bl.get_all_nodes()){
		auto nz = n->info().point.Z();
		ASSERT_TRUE((nz-5<1e-4) || (nz+5<1e-4) || (nz-2<1e-4));
	}
	ASSERT_EQ(18,bl.get_nb_cells<0>());
	ASSERT_EQ(33,bl.get_nb_cells<1>());
	ASSERT_EQ(20,bl.get_nb_cells<2>());
	ASSERT_EQ(4,bl.get_nb_cells<3>());
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, init_from_geom_bounding_box)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	ASSERT_EQ(8, bl.get_nb_cells<0>());
	ASSERT_EQ(12, bl.get_nb_cells<1>());
	ASSERT_EQ(6, bl.get_nb_cells<2>());
	ASSERT_EQ(1, bl.get_nb_cells<3>());

	for (auto a: bl.gmap()->attributes<0>()){
	       gmds::math::Point p = a.info().point;
			 ASSERT_NEAR(5, std::abs(p.X()),1e-8);
		    ASSERT_NEAR(5, std::abs(p.Y()),1e-8);
		    ASSERT_NEAR(5, std::abs(p.Z()),1e-8);
	}
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, single_block_to_mesh)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model);

	gmds::math::Point p000(0, 0, 0);
	gmds::math::Point p010(0, 1, 0);
	gmds::math::Point p110(1, 1, 0);
	gmds::math::Point p100(1, 0, 0);

	gmds::math::Point p001(0, 0, 1);
	gmds::math::Point p011(0, 1, 1);
	gmds::math::Point p111(1, 1, 1);
	gmds::math::Point p101(1, 0, 1);

	gmds::math::Point p002(0, 0, 2);
	gmds::math::Point p012(0, 1, 2);
	gmds::math::Point p112(1, 1, 2);
	gmds::math::Point p102(1, 0, 2);
	auto b = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);

	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
	bl.convert_to_mesh(m);

	ASSERT_EQ(m.getNbNodes(), 8);
	ASSERT_EQ(m.getNbEdges(), 12);
	ASSERT_EQ(m.getNbFaces(), 6);
	ASSERT_EQ(m.getNbRegions(), 1);
}


/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, projection_point_to_edges)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model, true);

	std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge> > all_edges = bl.get_all_sheet_edge_sets();
	gmds::math::Point p(5,5,5);
	for(auto sh_edges: all_edges){
		    ASSERT_EQ(4, sh_edges.size());
		    auto dist_coord = bl.get_projection_info(p,sh_edges);
		    double min_dist = MAXFLOAT;
		    for (auto dc: dist_coord){
				if (dc.first<min_dist)
				   min_dist = dc.first;
		    }
		    ASSERT_NEAR(min_dist,0,0.01);
	}
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_topological_queries)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model, true);
	std::vector<gmds::blocking::CurvedBlocking::Node> bl_nodes = bl.get_all_nodes();
	std::vector<gmds::blocking::CurvedBlocking::Edge> bl_edges = bl.get_all_edges();
	std::vector<gmds::blocking::CurvedBlocking::Face> bl_faces = bl.get_all_faces();

	for(auto n: bl_nodes){
		    std::vector<gmds::blocking::CurvedBlocking::Face> fs = bl.get_faces_of_node(n);
		    ASSERT_EQ(3, fs.size());
		    std::vector<gmds::blocking::CurvedBlocking::Edge> es = bl.get_edges_of_node(n);
		    ASSERT_EQ(3, es.size());
		    std::vector<gmds::blocking::CurvedBlocking::Block> bs = bl.get_blocks_of_node(n);
		    ASSERT_EQ(1, bs.size());
	}
	for(auto e: bl_edges){
		    std::vector<gmds::blocking::CurvedBlocking::Node> ns = bl.get_nodes_of_edge(e);
		    ASSERT_EQ(2, ns.size());
		    std::vector<gmds::blocking::CurvedBlocking::Face> fs = bl.get_faces_of_edge(e);
		    ASSERT_EQ(2, fs.size());
		    std::vector<gmds::blocking::CurvedBlocking::Block> bs = bl.get_blocks_of_edge(e);
		    ASSERT_EQ(1, bs.size());
	}
	for(auto f: bl_faces){
		    std::vector<gmds::blocking::CurvedBlocking::Node> ns = bl.get_nodes_of_face(f);
		    ASSERT_EQ(4, ns.size());
		    std::vector<gmds::blocking::CurvedBlocking::Edge> es = bl.get_edges_of_face(f);
		    ASSERT_EQ(4, es.size());
		    std::vector<gmds::blocking::CurvedBlocking::Block> bs = bl.get_blocks_of_face(f);
		    ASSERT_EQ(1, bs.size());
	}
}
TEST(CurvedBlockingTestSuite, save_vtk_blocking){
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);

	bl.save_vtk_blocking("testSaveWork.vtk");
}