#include <gmds/morphMesh/FastLocalize.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
/*                               TESTS UNITAIRES                              */
/*----------------------------------------------------------------------------*/

TEST(MorphMesh_FastLocalizeTestSuite, test_FastLocalize_2D)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::F | gmds::F2N));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/2D/C1_2D_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	// fastlocalize works on a point cloud, no need for the faces
   for(auto i: m.faces()) {
		m.deleteFace(i);
	}

	gmds::FastLocalize fl(&m);

	gmds::Cell::Data data = fl.find(gmds::math::Point({-0.0338,0.74,0}));
	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 734);

	data = fl.find(gmds::math::Point({-0.2,0.46,0}));
	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 11);
}


TEST(MorphMesh_FastLocalizeTestSuite, test_FastLocalize_3D)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::F | gmds::F2N));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/3D/C1_3D_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	// fastlocalize works on a point cloud, no need for the faces
	for(auto i: m.faces()) {
		m.deleteFace(i);
	}

	gmds::FastLocalize fl(&m);

	gmds::Cell::Data data = fl.find(gmds::math::Point({-1.0,0.0,0}));
	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 14322);

	data = fl.find(gmds::math::Point({-1.33,0.678,-1.89}));
	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 4048);
}


TEST(MorphMesh_FastLocalizeTestSuite, test_FastLocalize_nodes)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::F | gmds::F2N));

	for(int i=0; i<50; i++) {
		m.newNode(i,1,1);
	}

	gmds::FastLocalize fl(&m);

	gmds::Cell::Data data = fl.find(gmds::math::Point({1.9,0.0,0}));
	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 2);
}