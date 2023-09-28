//
// Created by rochec on 26/09/23.
//

#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/claire/MFEMMeshWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(MFEMMeshWriterTestClass, MFEM_Test)
{
	// Mesh
	Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));
	std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir+"/Aero/2D/APOLLO_2D_toFit.vtk";
	//std::string vtk_file = dir+"/Aero/2D/APOLLO_2D_20k_MESH.vtk";
	std::string vtk_file = dir+"/Aero/2D/Test2.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	//MFEMMeshWriter MeshWriter(&m, "APOLLO_2D_toFit");
	//MFEMMeshWriter MeshWriter(&m, "APOLLO_2D_20k_MESH");
	MFEMMeshWriter MeshWriter(&m, "Test2");
	MFEMMeshWriter::STATUS res = MeshWriter.execute();

	ASSERT_EQ(res, MeshWriter.SUCCESS);

}
