//
// Created by rochec on 19/05/2022.
//

#ifndef GMDS_SU2WRITERTESTSUITE_H
#define GMDS_SU2WRITERTESTSUITE_H

/*----------------------------------------------------------------------------*/
#include <gmds/aero/SU2Writer.h>
#include <gmds/aero/SU2Writer_3D.h>
#include <gmds/aero/AeroBoundaries_2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/aero/Utils.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(AeroTestClass, testSU2Writer)
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
	vtkReader.read(dir+"/Aero/2D/C1_2D_0.5.vtk");

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	math::Utils::MeshCleaner(&m);

	std::cout << "-> Ecriture du maillage ..." << std::endl;
	SU2Writer writer(&m, "SU2Writer_Test.su2", 1);
	SU2Writer::STATUS result = writer.execute();

	/*
	ioService = &m;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("Test.vtk");
	 */

	ASSERT_EQ(SU2Writer::SUCCESS, result);
}


/*
TEST(AeroTestClass, SU2Writer)
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
	vtkReader.read(dir+"/Aero/2D/Diamond_Airfoil_2D_Papier_5.vtk");

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	math::Utils::MeshCleaner(&m);

	std::cout << "-> Ecriture du maillage ..." << std::endl;
	SU2Writer writer(&m, "SU2Writer_Test.su2", -1000);
	SU2Writer::STATUS result = writer.execute();

	ASSERT_EQ(SU2Writer::SUCCESS, result);
}
*/


#endif     // GMDS_SU2WRITERTESTSUITE_H
