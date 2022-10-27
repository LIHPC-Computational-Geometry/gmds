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

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReader.read(dir + "/Aero/2D/Apollo_2D_5.vtk");

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();

	math::Utils::MeshCleaner(&m);

	AeroBoundaries_2D* Bnd = new AeroBoundaries_2D(&m) ;
	Bnd->execute();

	int mark_Farfiel = Bnd->getMarkAmont();
	int mark_Paroi = Bnd->getMarkParoi();

	Variable<double>* var_ls_test = m.newVariable<double,GMDS_NODE>("GMDS_Distance_TEST");
	DiffusionEquation2D ls_test(&m, mark_Paroi, mark_Farfiel, var_ls_test);
	DiffusionEquation2D::STATUS ls_result = ls_test.execute();

	ASSERT_EQ(DiffusionEquation2D::SUCCESS, ls_result);

	double eps = pow(10, -5);

	{
		ASSERT_NEAR( var_ls_test->value(0), -1.14684e-15, eps );
		ASSERT_NEAR( var_ls_test->value(100), 1, eps );
		ASSERT_NEAR( var_ls_test->value(200), 1, eps );
		ASSERT_NEAR( var_ls_test->value(300), 0.0635951, eps );
		ASSERT_NEAR( var_ls_test->value(400), 0.499526, eps );
		ASSERT_NEAR( var_ls_test->value(500), 0.635191, eps );
		ASSERT_NEAR( var_ls_test->value(600), 0.645018, eps );
		ASSERT_NEAR( var_ls_test->value(700), 0.474929, eps );
		ASSERT_NEAR( var_ls_test->value(800), 0.660854, eps );
		ASSERT_NEAR( var_ls_test->value(900), 0.752721, eps );
		ASSERT_NEAR( var_ls_test->value(1000), 0.79844, eps );
		ASSERT_NEAR( var_ls_test->value(1100), 0.329108, eps );
		ASSERT_NEAR( var_ls_test->value(1200), 0.859997, eps );
		ASSERT_NEAR( var_ls_test->value(1300), 0.853981, eps );
		ASSERT_NEAR( var_ls_test->value(1400), 0.86001, eps );
		ASSERT_NEAR( var_ls_test->value(1500), 0.444049, eps );
		ASSERT_NEAR( var_ls_test->value(1600), 0.734511, eps );
		ASSERT_NEAR( var_ls_test->value(1700), 0.894059, eps );
		ASSERT_NEAR( var_ls_test->value(1800), 0.155853, eps );
		ASSERT_NEAR( var_ls_test->value(1900), 0.108677, eps );
		ASSERT_NEAR( var_ls_test->value(2000), 0.536911, eps );
		ASSERT_NEAR( var_ls_test->value(2100), 0.76033, eps );
		ASSERT_NEAR( var_ls_test->value(2200), 0.69826, eps );
		ASSERT_NEAR( var_ls_test->value(2300), 0.501803, eps );
		ASSERT_NEAR( var_ls_test->value(2400), 0.845343, eps );
		ASSERT_NEAR( var_ls_test->value(2500), 0.643272, eps );
		ASSERT_NEAR( var_ls_test->value(2600), 0.901117, eps );
		ASSERT_NEAR( var_ls_test->value(2700), 0.939982, eps );
		ASSERT_NEAR( var_ls_test->value(2800), 0.969981, eps );
		ASSERT_NEAR( var_ls_test->value(2900), 0.852891, eps );
		ASSERT_NEAR( var_ls_test->value(3000), 0.119069, eps );
		ASSERT_NEAR( var_ls_test->value(3100), 0.901179, eps );
	}

}
