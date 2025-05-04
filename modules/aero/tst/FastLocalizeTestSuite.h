//
// Created by rochec on 23/03/2022.
//

#include <gmds/aero/FastLocalize.h>
#include <gmds/ig/Mesh.h>
#include <gmds/aero/Utils.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*                               TESTS UNITAIRES                              */
/*----------------------------------------------------------------------------*/

TEST(FastLocalizeTestSuite, test_FastLocalize_2D)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/2D/C1_2D_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	FastLocalize fl(&m);
	gmds::Cell::Data data = fl.find(gmds::math::Point({-0.0338,0.74,0}));

	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 734);
	data = fl.find(gmds::math::Point({-0.2,0.46,0}));

	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 11);

}


TEST(FastLocalizeTestSuite, test_FastLocalize_3D)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/3D/C1_3D_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	FastLocalize fl(&m);

	gmds::Cell::Data data = fl.find(gmds::math::Point({-1.0,0.0,0}));
	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 14322);

	data = fl.find(gmds::math::Point({-1.33,0.678,-1.89}));
	ASSERT_EQ(data.dim, 0);
	ASSERT_EQ(data.id, 4048);

}