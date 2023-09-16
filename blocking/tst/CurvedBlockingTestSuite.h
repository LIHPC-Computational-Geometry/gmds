/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
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
std::tuple<int,int,int,int> get_node_statistics(gmds::blocking::CurvedBlocking& ABlocking){
    auto nb_on_vertex=0;
    auto nb_on_curve=0;
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gmds::blocking::CurvedBlocking::Node> all_nodes = ABlocking.get_all_nodes();
    for(auto n:all_nodes){
        if(n->info().geom_dim==0)
            nb_on_vertex++;
        else if(n->info().geom_dim==1)
            nb_on_curve++;
        else if(n->info().geom_dim==2)
            nb_on_surface++;
        else if(n->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_vertex,nb_on_curve,nb_on_surface,nb_in_volume);
}
/*----------------------------------------------------------------------------*/
std::tuple<int,int,int> get_edge_statistics(gmds::blocking::CurvedBlocking& ABlocking){
    auto nb_on_curve=0;
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gmds::blocking::CurvedBlocking::Edge> all_edges = ABlocking.get_all_edges();
    for(auto e:all_edges){
        if(e->info().geom_dim==1)
            nb_on_curve++;
        else if(e->info().geom_dim==2)
            nb_on_surface++;
        else if(e->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_curve,nb_on_surface,nb_in_volume);
}
/*----------------------------------------------------------------------------*/
std::tuple<int,int> get_face_statistics(gmds::blocking::CurvedBlocking& ABlocking){
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = ABlocking.get_all_faces();
    for(auto f:all_faces){
        if(f->info().geom_dim==2)
            nb_on_surface++;
        else if(f->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_surface,nb_in_volume);
}
/*----------------------------------------------------------------------------*/
void export_vtk(gmds::blocking::CurvedBlocking& ABlocking, int AModel, const std::string& AFileName){
    gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
    ABlocking.convert_to_mesh(m_out);
    gmds::IGMeshIOService ioService(&m_out);
    gmds::VTKWriter writer(&ioService);
    writer.setCellOptions(AModel);
    writer.setDataOptions(AModel);
    writer.write(AFileName);
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
	ASSERT_EQ(b->info().geom_dim, 3);
	ASSERT_EQ(b->info().geom_id, gmds::NullID);
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

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_init_from_ig_mesh)
{
	gmds::cad::FACManager geom_model;
	setUp(geom_model);
	gmds::blocking::CurvedBlocking bl(&geom_model, false);

	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N |  gmds::R | gmds::R2N));
	gmds::Node n0 = m.newNode(gmds::math::Point(0,0,0));
	gmds::Node n1 = m.newNode(gmds::math::Point(1,0,0));
	gmds::Node n2 = m.newNode(gmds::math::Point(1,1,0));
	gmds::Node n3 = m.newNode(gmds::math::Point(0,1,0));

	gmds::Node n4 = m.newNode(gmds::math::Point(0,0,1));
	gmds::Node n5 = m.newNode(gmds::math::Point(1,0,1));
	gmds::Node n6 = m.newNode(gmds::math::Point(1,1,1));
	gmds::Node n7 = m.newNode(gmds::math::Point(0,1,1));

	gmds::Node n8 = m.newNode(gmds::math::Point(0,0,2));
	gmds::Node n9 = m.newNode(gmds::math::Point(1,0,2));
	gmds::Node n10= m.newNode(gmds::math::Point(1,1,2));
	gmds::Node n11= m.newNode(gmds::math::Point(0,1,2));


	gmds::Node n12= m.newNode(gmds::math::Point(0,0,3));
	gmds::Node n13= m.newNode(gmds::math::Point(1,0,3));
	gmds::Node n14= m.newNode(gmds::math::Point(1,1,3));
	gmds::Node n15= m.newNode(gmds::math::Point(0,1,3));

	m.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
	m.newHex(n4,n5,n6,n7,n8,n9,n10,n11);
	m.newHex(n8,n9,n10,n11,n12,n13,n14,n15);
	bl.init_from_mesh(m);

	ASSERT_EQ(16,bl.get_all_nodes().size());
	ASSERT_EQ(28,bl.get_all_edges().size());
	ASSERT_EQ(16,bl.get_all_faces().size());
	ASSERT_EQ(3,bl.get_all_blocks().size());
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_chord_query)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, false);

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N |  gmds::R | gmds::R2N));
    gmds::Node n0 = m.newNode(gmds::math::Point(0,0,0));
    gmds::Node n1 = m.newNode(gmds::math::Point(1,0,0));
    gmds::Node n2 = m.newNode(gmds::math::Point(1,1,0));
    gmds::Node n3 = m.newNode(gmds::math::Point(0,1,0));

    gmds::Node n4 = m.newNode(gmds::math::Point(0,0,1));
    gmds::Node n5 = m.newNode(gmds::math::Point(1,0,1));
    gmds::Node n6 = m.newNode(gmds::math::Point(1,1,1));
    gmds::Node n7 = m.newNode(gmds::math::Point(0,1,1));

    gmds::Node n8 = m.newNode(gmds::math::Point(0,0,2));
    gmds::Node n9 = m.newNode(gmds::math::Point(1,0,2));
    gmds::Node n10= m.newNode(gmds::math::Point(1,1,2));
    gmds::Node n11= m.newNode(gmds::math::Point(0,1,2));

    gmds::Node n12= m.newNode(gmds::math::Point(0,0,3));
    gmds::Node n13= m.newNode(gmds::math::Point(1,0,3));
    gmds::Node n14= m.newNode(gmds::math::Point(1,1,3));
    gmds::Node n15= m.newNode(gmds::math::Point(0,1,3));

    m.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
    m.newHex(n4,n5,n6,n7,n8,n9,n10,n11);
    m.newHex(n8,n9,n10,n11,n12,n13,n14,n15);
    bl.init_from_mesh(m);


    std::vector<gmds::blocking::Dart3> darts;
    std::vector<gmds::blocking::CurvedBlocking::Face> bl_faces = bl.get_all_faces();
    for(auto f:bl_faces){
        std::cout<<bl.get_center_of_face(f)<<std::endl;
        gmds::math::Point center = bl.get_center_of_face(f);
        double z_integer_part = std::floor(center.Z());
        double z_decimal_part = center.Z()-z_integer_part;
        if(z_decimal_part==0){
            bl.get_all_chord_darts(f, darts);
            ASSERT_EQ(darts.size(),4);
            ASSERT_EQ(bl.get_all_chord_blocks(f).size(),3);
        }
        else{
            bl.get_all_chord_darts(f, darts);
            ASSERT_EQ(darts.size(),2);
            ASSERT_EQ(bl.get_all_chord_blocks(f).size(),1);
        }
    }
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_chord_collapse)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, false);

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

    gmds::GridBuilder gb(&m,3);
    gb.execute(4,1.0, 4, 1.0, 4, 1.0);
    bl.init_from_mesh(m);


    std::vector<gmds::blocking::CurvedBlocking::Face> bl_faces = bl.get_all_faces();

    ASSERT_EQ(bl.get_nb_cells<3>(),27);

    bool found_face = false;
    auto face_id = -1;
    gmds::math::Point seed(1.5,1.5,0.0);
    for(auto i=0; i<bl_faces.size() && !found_face;i++){
        gmds::blocking::CurvedBlocking::Face fi = bl_faces[i];
        gmds::math::Point ci = bl.get_center_of_face(fi);
        if(ci.distance2(seed)<0.1){
            found_face = true;
            face_id=i;
        }
    }

    std::vector<gmds::blocking::CurvedBlocking::Node> f_nodes = bl.get_nodes_of_face(bl_faces[face_id]);
    bl.collapse_chord(bl_faces[face_id],f_nodes[0],f_nodes[2]);

    ASSERT_EQ(bl.get_nb_cells<3>(),24);

    gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
    bl.convert_to_mesh(m_out);

    gmds::IGMeshIOService ioService(&m_out);
    gmds::VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N | gmds::F);
    writer.setDataOptions(gmds::N |gmds::F );
    writer.write("collapse.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_validate_pillow)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    gmds::cad::GeomVolume* v = geom_model.getVolumes()[0];

    auto [x_min,y_min,z_min,x_max,y_max,z_max] = v->BBox();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);
    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;
    //We pick all the boundary faces
    for(auto f:all_faces){
        if(f->info().geom_dim==2){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.validate_pillowing_surface(surf));

    surf.clear();
    //We pick a full boundary surface on coord X=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.validate_pillowing_surface(surf));
    surf.pop_back();
    ASSERT_FALSE(bl.validate_pillowing_surface(surf));
    surf.pop_back();
    ASSERT_FALSE(bl.validate_pillowing_surface(surf));
    surf.pop_back();
    ASSERT_FALSE(bl.validate_pillowing_surface(surf));
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    ASSERT_EQ(bl.get_nb_cells<3>(),1);
    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(1);
//    export_vtk(bl,gmds::N | gmds::F, "pillow_1_surf.vtk");
    ASSERT_EQ(bl.get_nb_cells<3>(),2);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 4);
    ASSERT_EQ(nb_nodes_on_surface, 0);
    ASSERT_EQ(nb_nodes_in_volume, 0);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 16);
    ASSERT_EQ(nb_edges_on_surface, 4);
    ASSERT_EQ(nb_edges_in_volume, 0);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 10);
    ASSERT_EQ(nb_faces_in_volume, 1);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_2)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }
    ASSERT_EQ(bl.get_nb_cells<3>(),8);
    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;
    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_2_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),16);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 16);
    ASSERT_EQ(nb_nodes_on_surface, 14);
    ASSERT_EQ(nb_nodes_in_volume, 4);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 28);
    ASSERT_EQ(nb_edges_on_surface, 44);
    ASSERT_EQ(nb_edges_in_volume, 19);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 36);
    ASSERT_EQ(nb_faces_in_volume, 30);
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_3)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_3_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),3);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 4);
    ASSERT_EQ(nb_nodes_on_surface, 2);
    ASSERT_EQ(nb_nodes_in_volume, 0);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 16);
    ASSERT_EQ(nb_edges_on_surface, 8);
    ASSERT_EQ(nb_edges_in_volume, 1);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 12);
    ASSERT_EQ(nb_faces_in_volume, 3);}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_4)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);
    //export_vtk(bl,gmds::N | gmds::F, "pillow_4_surf.vtk");
    ASSERT_EQ(bl.get_nb_cells<3>(),4);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 3);
    ASSERT_EQ(nb_nodes_on_surface, 3);
    ASSERT_EQ(nb_nodes_in_volume, 1);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 15);
    ASSERT_EQ(nb_edges_on_surface, 9);
    ASSERT_EQ(nb_edges_in_volume, 4);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 12);
    ASSERT_EQ(nb_faces_in_volume, 6);
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_5)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()+5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);
    //export_vtk(bl,gmds::N | gmds::F, "pillow_5_surf.vtk");
    ASSERT_EQ(bl.get_nb_cells<3>(),5);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 2);
    ASSERT_EQ(nb_nodes_on_surface, 4);
    ASSERT_EQ(nb_nodes_in_volume, 2);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 14);
    ASSERT_EQ(nb_edges_on_surface, 10);
    ASSERT_EQ(nb_edges_in_volume, 7);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 12);
    ASSERT_EQ(nb_faces_in_volume, 9);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_6)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()+5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()+5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_6_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),6);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 0);
    ASSERT_EQ(nb_nodes_on_surface, 4);
    ASSERT_EQ(nb_nodes_in_volume, 4);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 12);
    ASSERT_EQ(nb_edges_on_surface, 8);
    ASSERT_EQ(nb_edges_in_volume, 12);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 10);
    ASSERT_EQ(nb_faces_in_volume, 13);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_7)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();

    ASSERT_TRUE(bl.pillow(all_faces));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_7_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),7);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 0);
    ASSERT_EQ(nb_nodes_on_surface, 0);
    ASSERT_EQ(nb_nodes_in_volume, 8);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 12);
    ASSERT_EQ(nb_edges_on_surface, 0);
    ASSERT_EQ(nb_edges_in_volume, 20);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 6);
    ASSERT_EQ(nb_faces_in_volume, 18);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_8)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0
    for(auto f:all_faces){
        auto nb_adj = bl.get_blocks_of_face(f);
        if(nb_adj.size()==1)
            surf.push_back(f);
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_8_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),32);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 12);
    ASSERT_EQ(nb_nodes_on_surface, 6);
    ASSERT_EQ(nb_nodes_in_volume, 27);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 24);
    ASSERT_EQ(nb_edges_on_surface, 24);
    ASSERT_EQ(nb_edges_in_volume, 80);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 24);
    ASSERT_EQ(nb_faces_in_volume, 84);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_9)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_9_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),20);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 15);
    ASSERT_EQ(nb_nodes_on_surface, 15);
    ASSERT_EQ(nb_nodes_in_volume, 8);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 27);
    ASSERT_EQ(nb_edges_on_surface, 45);
    ASSERT_EQ(nb_edges_in_volume, 31);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 36);
    ASSERT_EQ(nb_faces_in_volume, 42);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_10)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick an inner surface
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X())<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_EQ(surf.size(),4);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_10_surf_before.vtk");
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    export_vtk(bl,gmds::N | gmds::E, "pillow_10_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),12);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 16);
    ASSERT_EQ(nb_nodes_on_surface, 10);
    ASSERT_EQ(nb_nodes_in_volume, 2);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 28);
    ASSERT_EQ(nb_edges_on_surface, 36);
    ASSERT_EQ(nb_edges_in_volume, 11);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 32);
    ASSERT_EQ(nb_faces_in_volume, 20);
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_11)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick an inner surface
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X())<0.1 && ci.Y()<0 && ci.Z()<0){
            surf.push_back(f);
        }
        else if(ci.X() <0 && fabs(ci.Y())<0.1 && ci.Z()<0){
            surf.push_back(f);
        }
        else if(ci.X()<0 && ci.Y()<0 && fabs(ci.Z())<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_EQ(surf.size(),3);

    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

   // export_vtk(bl,gmds::N | gmds::E, "pillow_11_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),11);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 15);
    ASSERT_EQ(nb_nodes_on_surface, 9);
    ASSERT_EQ(nb_nodes_in_volume, 2);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 27);
    ASSERT_EQ(nb_edges_on_surface, 33);
    ASSERT_EQ(nb_edges_in_volume, 10);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 30);
    ASSERT_EQ(nb_faces_in_volume, 18);
}

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, test_pillow_12)
{
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick an inner surface
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X())<0.1 ){
            surf.push_back(f);
        }
        else if(ci.X() <0 && fabs(ci.Y()-5)<0.1 ){
            surf.push_back(f);
        }
    }
    ASSERT_EQ(surf.size(),6);

    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::E, "pillow_12_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),14);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 16);
    ASSERT_EQ(nb_nodes_on_surface, 12);
    ASSERT_EQ(nb_nodes_in_volume, 3);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 28);
    ASSERT_EQ(nb_edges_on_surface, 40);
    ASSERT_EQ(nb_edges_in_volume, 15);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 34);
    ASSERT_EQ(nb_faces_in_volume, 25);
}