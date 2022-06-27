//
// Created by rochec on 10/02/2022.
//

#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/AeroPipeline_2D.h>
#include <gmds/claire/AeroPipeline_3D.h>
#include <gmds/claire/AeroExtrusion_2D.h>
#include <gmds/claire/AeroException.h>
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

	// Paramètres de dimension
	params_aero.dim=ParamsAero::DIM_2D;

	// Paramètres d'entrée/sortie
	std::string dir(TEST_SAMPLES_DIR);
	params_aero.input_file=dir+"/Aero/Poubelle/mesh_1.vtk";
	params_aero.output_file="AeroPipeline2D_Quad.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	// Paramètres discrétisation de la paroi
	params_aero.nbrMinBloc=16;

	// Paramètres de l'extrusion
	params_aero.vectors_field = 0;
	params_aero.nbr_couches = 10;
	params_aero.delta_cl = 0.025;			// Epaisseur de la première couche
	params_aero.x_lim=0;			// Limites physiques à partir desquelles
	params_aero.y_lim=-10000;			// l'insertion et la fusion de blocs
	params_aero.z_lim=-10000;			// sont autorisées

	AeroPipeline_2D algo_aero2D(params_aero);
	AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero2D_result);

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
	params_aero.input_file=dir+"/Aero/3D/biconique.vtk";
	params_aero.output_file="AeroPipeline3D_Hexa.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	AeroPipeline_3D algo_aero3D(params_aero);
	AbstractAeroPipeline::STATUS aero3D_result = algo_aero3D.execute();

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero3D_result);

}


/*----------------------------------------------------------------------------*/
/*              CAS TEST 3D CLASSE LevelSet avec Exceptions                   */
/*----------------------------------------------------------------------------*/

/*
TEST(AeroPipelineTestClass, AeroExtrusion_Test1)
{
	AeroExtrusion_2D extrusion2D(NULL);
	try {
		extrusion2D.execute();
	}
	catch (AeroException& e) {

		std::cout<<"Exception catched"<<std::endl;
	}
	catch (...){
		std::cout<<"any exception"<<std::endl;
	}
	std::cout<<"here"<<std::endl;
}
*/