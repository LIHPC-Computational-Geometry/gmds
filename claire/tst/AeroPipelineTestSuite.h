//
// Created by rochec on 10/02/2022.
//

#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/AeroPipeline2D.h>
#include <gmds/claire/AeroPipeline3D.h>
#include <gmds/claire/Params.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*                       CAS TEST 2D CLASSE LevelSet                          */
/*----------------------------------------------------------------------------*/

TEST(AeroPipelineTestClass, AeroPipeline2D_Test1)
{
	ParamsAero params_aero;

	// Définition en dur des paramètres pour l'aéro
	// le temps de mettre en place un fichier .ini
	params_aero.dim=ParamsAero::DIM_2D;
	std::string dir(TEST_SAMPLES_DIR);
	params_aero.input_file=dir+"/Aero/1a_apollo_experimental_coupe_5.vtk";
	params_aero.output_file="AeroPipeline2D_Quad.vtk";
	params_aero.output_dir="gmds/claire/tst/";
	params_aero.nbrMinBloc=16;

	AeroPipeline2D algo_aero2D(params_aero);
	algo_aero2D.execute();

	bool isOver = algo_aero2D.getIsOver();

	ASSERT_EQ(isOver, true);

}



/*----------------------------------------------------------------------------*/
/*                       CAS TEST 3D CLASSE LevelSet                          */
/*----------------------------------------------------------------------------*/

TEST(AeroPipelineTestClass, AeroPipeline3D_Test1)
{
	ParamsAero params_aero;

	// Définition en dur des paramètres pour l'aéro
	// le temps de mettre en place un fichier .ini
	params_aero.dim=ParamsAero::DIM_3D;
	std::string dir(TEST_SAMPLES_DIR);
	params_aero.input_file=dir+"/Aero/3D/caretwing.vtk";
	params_aero.output_file="AeroPipeline3D_Test1_Result.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	AeroPipeline3D algo_aero3D(params_aero);
	algo_aero3D.execute();

	bool isOver = algo_aero3D.getIsOver();

	ASSERT_EQ(isOver, true);

}
