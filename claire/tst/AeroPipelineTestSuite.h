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

	// Dimension Parameter
	params_aero.dim=ParamsAero::DIM_2D;

	// IN/OUT Parameters
	std::string dir(TEST_SAMPLES_DIR);
	params_aero.input_file=dir+"/Aero/2D/NACA_0012_2D_0.5.vtk";
	params_aero.output_file="AeroPipeline2D_Quad.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	// Physical Parameters for the algorithm
	params_aero.delta_cl = 0.08;				// Epaisseur de la première couche, pour la couche limite
	params_aero.angle_attack = 0;			// Angle of attack (in degrees)

	// Wall Discretization Parameter
	params_aero.nbrMinBloc = 10;					// Minimal number of block on the wall
	params_aero.nbrCellsInCL = 30;				// Number of cells in the boundary layer
	params_aero.cell_size_dx_wall = 0.1 ;			// Maximal size of the cells on the wall
	params_aero.cell_size_default = 0.08 ;			// Default size of the cells in the domain

	// Vector Field Computation Parameters
	params_aero.vectors_field = 0;				// Choose the way the vectors field is computed for the extrusion
	params_aero.x_VectorField_Z1 = 0.4;			// Choose the x value of the first zone  [-inf, x_VectorField_Z1]
	params_aero.x_VectorField_Z2 = 1.4;			// Choose the x value of the second zone [x_VectorField_Z2, +inf]

	// Extrusion Parameters
	params_aero.nbr_couches = 3;			// Number of layer in extrusion
	params_aero.x_lim = 1.0;			// Limites physiques à partir desquelles
	params_aero.y_lim = -10000;			// l'insertion et la fusion de blocs
	params_aero.z_lim = -10000;			// sont autorisées


	// SU2 Writer Parameter
	params_aero.x_lim_SU2_inoutlet = -pow(10,6);		// Limit between inlet and outlet for SU2 writer

	// Mesh Generation
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