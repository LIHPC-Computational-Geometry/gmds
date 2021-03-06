//
// Created by rochec on 19/05/2022.
//

#ifndef GMDS_SU2WRITERTESTSUITE_H
#define GMDS_SU2WRITERTESTSUITE_H

/*----------------------------------------------------------------------------*/
#include <gmds/claire/SU2Writer.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
/*TEST(ClaireTestClass, testSU2Writer)
{
	// Test
	Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                       gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	std::string dir(TEST_SAMPLES_DIR);
	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(dir+"/Aero/2D/SU2_C1_2D_0.1.vtk");

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	SU2Writer writer(&m, "Test.su2", -10.0);
	SU2Writer::STATUS result = writer.execute();

	ioService = &m;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("Test.vtk");

	ASSERT_EQ(SU2Writer::SUCCESS, result);
}
*/
#endif     // GMDS_SU2WRITERTESTSUITE_H
