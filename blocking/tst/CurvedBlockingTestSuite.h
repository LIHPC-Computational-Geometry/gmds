/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
void
setUp(gmds::cad::FACManager &AGeomManager)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
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
TEST(CurvedBlockingTestSuite, init)
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
	ASSERT_EQ(2, bl.get_nb_cells<3>());
	bl.sew<3>(b1->dart(), b2->dart());
	ASSERT_EQ(12, bl.get_nb_cells<0>());
	ASSERT_EQ(20, bl.get_nb_cells<1>());
	ASSERT_EQ(11, bl.get_nb_cells<2>());
	ASSERT_EQ(2, bl.get_nb_cells<3>());
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