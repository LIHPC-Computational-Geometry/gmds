//
// Created by rochec on 25/10/2022.
//

#include <gmds/claire/Utils.h>
#include <gmds/claire/DiffusionEquation2D.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(DiffusionEquation2DTestClass, DiffusionEquation2D)
{

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/2D/Apollo_2D_5.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	math::Utils::MeshCleaner(&m);


	AeroBoundaries_2D* Bnd = new AeroBoundaries_2D(&m) ;
	Bnd->execute();

	int mark_Farfiel = Bnd->getMarkAmont();
	int mark_Paroi = Bnd->getMarkParoi();

	Variable<double>* var_ls_test = m.newVariable<double,GMDS_NODE>("GMDS_Distance_TEST");
	DiffusionEquation2D ls_test(&m, mark_Paroi, mark_Farfiel, var_ls_test);
	DiffusionEquation2D::STATUS ls_result = ls_test.execute();

	ASSERT_EQ(DiffusionEquation2D::SUCCESS, ls_result);

}
