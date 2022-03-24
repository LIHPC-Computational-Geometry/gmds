//
// Created by rochec on 23/03/2022.
//

#include <gmds/claire/AbstractAeroBoundaries.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/io/IGMeshIOService.h>
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

TEST(AeroBoundariesTestClass, AeroBoundaries2D_Test1)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/C1_2D_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	// Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	AeroBoundaries_2D bnd_2D(&m);
	AbstractAeroBoundaries::STATUS bnd_2D_result = bnd_2D.execute();

	ASSERT_EQ(bnd_2D_result, AbstractAeroBoundaries::SUCCESS);

	ASSERT_EQ(bnd_2D.isParoi(6), true);
	ASSERT_EQ(bnd_2D.isParoi(7), true);
	ASSERT_EQ(bnd_2D.isParoi(8), true);
	ASSERT_EQ(bnd_2D.isParoi(9), true);
	ASSERT_EQ(bnd_2D.isParoi(10), true);
	ASSERT_EQ(bnd_2D.isParoi(3), true);
	ASSERT_EQ(bnd_2D.isParoi(32), true);
	ASSERT_EQ(bnd_2D.isParoi(1281), false);
	ASSERT_EQ(bnd_2D.isParoi(850), false);
	ASSERT_EQ(bnd_2D.isParoi(872), false);

	ASSERT_EQ(bnd_2D.isAmont(40), true);
	ASSERT_EQ(bnd_2D.isAmont(43), true);
	ASSERT_EQ(bnd_2D.isAmont(45), true);
	ASSERT_EQ(bnd_2D.isAmont(108), true);
	ASSERT_EQ(bnd_2D.isAmont(111), true);
	ASSERT_EQ(bnd_2D.isAmont(105), true);
	ASSERT_EQ(bnd_2D.isAmont(10), false);
	ASSERT_EQ(bnd_2D.isAmont(3), false);
	ASSERT_EQ(bnd_2D.isAmont(32), false);
	ASSERT_EQ(bnd_2D.isAmont(1018), false);
	ASSERT_EQ(bnd_2D.isAmont(550), false);

	ASSERT_EQ(bnd_2D.isBnd(40), true);
	ASSERT_EQ(bnd_2D.isBnd(43), true);
	ASSERT_EQ(bnd_2D.isBnd(111), true);
	ASSERT_EQ(bnd_2D.isBnd(29), true);
	ASSERT_EQ(bnd_2D.isBnd(32), true);
	ASSERT_EQ(bnd_2D.isBnd(13), true);
	ASSERT_EQ(bnd_2D.isBnd(16), true);
	ASSERT_EQ(bnd_2D.isBnd(17), true);
	ASSERT_EQ(bnd_2D.isBnd(1502), false);
	ASSERT_EQ(bnd_2D.isBnd(551), false);
	ASSERT_EQ(bnd_2D.isBnd(888), false);
	ASSERT_EQ(bnd_2D.isBnd(1361), false);
	ASSERT_EQ(bnd_2D.isBnd(278), false);
	ASSERT_EQ(bnd_2D.isBnd(400), false);
	ASSERT_EQ(bnd_2D.isBnd(744), false);
	ASSERT_EQ(bnd_2D.isBnd(418), false);

	ASSERT_EQ(bnd_2D.isImmerged(), true);
	ASSERT_EQ(bnd_2D.getNbrBords(), 2);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("AeroBoundaries_2D_Test1.vtk");


}
