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
	params_aero.input_file=dir+"/Aero/2D/Apollo_2D_5.vtk";
	params_aero.output_file="AeroPipeline2D_Quad.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	// Physical Parameters for the algorithm
	params_aero.delta_cl = 3;				// Epaisseur de la première couche, pour la couche limite
	params_aero.angle_attack = 0;			// Angle of attack (in degrees)

	// Wall Discretization Parameter
	params_aero.nbrMinBloc = 4;					// Minimal number of block on the wall
	params_aero.nbrCellsInCL = 30;				// Number of cells in the boundary layer
	params_aero.edge_size_wall = 3 ;										// Maximal size of the cells on the wall
	params_aero.edge_size_default = 5 ;										// Default size of the cells in the domain
	params_aero.edge_size_first_ortho_wall = 1*pow(10,-5);	// Size of the first edge orthogonal to the wall

	// Vector Field Computation Parameters
	params_aero.vectors_field = 1;				// Choose the way the vectors field is computed for the extrusion
	params_aero.x_VectorField_Z1 = 200;			// Choose the x value of the first zone  [-inf, x_VectorField_Z1]
	params_aero.x_VectorField_Z2 = 600;			// Choose the x value of the second zone [x_VectorField_Z2, +inf]

	// Extrusion Parameters
	params_aero.nbr_couches = 3;			// Number of layer in extrusion
	params_aero.x_lim = 100;				// Limites physiques à partir desquelles
	params_aero.y_lim = -10000;			// l'insertion et la fusion de blocs
	params_aero.z_lim = -10000;			// sont autorisées

	// Smoothing Parameters
	params_aero.nbr_iter_smoothing_yao = 0 ;		// Number of iterations for the Yao Smoothing
	params_aero.damping_smoothing_yao = 0.2 ;			// Damping parameter for the Yao Smoothing

	// SU2 Writer Parameter
	params_aero.x_lim_SU2_inoutlet = -pow(10,6);		// Limit between inlet and outlet for SU2 writer

	// Mesh Generation
	AeroPipeline_2D algo_aero2D(params_aero);
	AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero2D_result);

}

/*
TEST(AeroPipelineTestClass, AeroPipeline2D_User_NACA)
{
	ParamsAero params_aero;

	// Définition en dur des paramètres pour l'aéro
	// le temps de mettre en place un fichier .ini

	// Dimension Parameter
	params_aero.dim=ParamsAero::DIM_2D;

	// IN/OUT Parameters
	std::string dir(TEST_SAMPLES_DIR);
	params_aero.input_file=dir+"/Aero/2D/NACA_0012_2D_700_0.05.vtk";
	//params_aero.input_file=dir+"/Aero/2D/test_2.vtk";
	params_aero.output_file="AeroPipeline2D_Quad.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	// Physical Parameters for the algorithm
	params_aero.delta_cl = 0.04;				// Epaisseur de la première couche, pour la couche limite
	params_aero.angle_attack = 0;			// Angle of attack (in degrees)

	// Wall Discretization Parameter
	params_aero.nbrMinBloc = 33;					// Minimal number of block on the wall
	params_aero.nbrCellsInCL = 100;				// Number of cells in the boundary layer
	params_aero.edge_size_wall = 0.001 ;										// Maximal size of the cells on the wall
	params_aero.edge_size_default = 0.012 ;										// Default size of the cells in the domain
	params_aero.edge_size_first_ortho_wall = 1*pow(10,-8);	// Size of the first edge orthogonal to the wall

	// Vector Field Computation Parameters
	params_aero.vectors_field = 3;				// Choose the way the vectors field is computed for the extrusion
	params_aero.x_VectorField_Z1 = 1.5;			// Choose the x value of the first zone  [-inf, x_VectorField_Z1]
	params_aero.x_VectorField_Z2 = 5.0;			// Choose the x value of the second zone [x_VectorField_Z2, +inf]

	// Extrusion Parameters
	params_aero.nbr_couches = 4;			// Number of layer in extrusion
	params_aero.x_lim = 1;				// Limites physiques à partir desquelles
	params_aero.y_lim = -10000;			// l'insertion et la fusion de blocs
	params_aero.z_lim = -10000;			// sont autorisées

	// Smoothing Parameters
	params_aero.nbr_iter_smoothing_yao = 0 ;		// Number of iterations for the Yao Smoothing
	params_aero.damping_smoothing_yao = 0.2 ;		// Damping parameter for the Yao Smoothing

	// SU2 Writer Parameter
	params_aero.x_lim_SU2_inoutlet = -pow(10,6);		// Limit between inlet and outlet for SU2 writer

	// Mesh Generation
	AeroPipeline_2D algo_aero2D(params_aero);
	AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero2D_result);

}
*/

/*
TEST(AeroPipelineTestClass, AeroPipeline2D_User_DA)
{
	ParamsAero params_aero;

	// Définition en dur des paramètres pour l'aéro
	// le temps de mettre en place un fichier .ini

	// Dimension Parameter
	params_aero.dim=ParamsAero::DIM_2D;

	// IN/OUT Parameters
	std::string dir(TEST_SAMPLES_DIR);
	params_aero.input_file=dir+"/Aero/2D/Diamond_Airfoil_2D_Papier_150_5.vtk";
	//params_aero.input_file=dir+"/Aero/2D/test_2.vtk";
	params_aero.output_file="AeroPipeline2D_Quad.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	// Physical Parameters for the algorithm
	params_aero.delta_cl = 5.0;				// Epaisseur de la première couche, pour la couche limite
	params_aero.angle_attack = 0;			// Angle of attack (in degrees)

	// Wall Discretization Parameter
	params_aero.nbrMinBloc = 4;					// Minimal number of block on the wall
	params_aero.nbrCellsInCL = 30;				// Number of cells in the boundary layer
	params_aero.edge_size_wall = 0.4 ;										// Maximal size of the cells on the wall
	params_aero.edge_size_default = 1.0 ;										// Default size of the cells in the domain
	params_aero.edge_size_first_ortho_wall = 1*pow(10,-7);	// Size of the first edge orthogonal to the wall

	// Vector Field Computation Parameters
	params_aero.vectors_field = 3;				// Choose the way the vectors field is computed for the extrusion
	params_aero.x_VectorField_Z1 = 150;			// Choose the x value of the first zone  [-inf, x_VectorField_Z1]
	params_aero.x_VectorField_Z2 = 600;			// Choose the x value of the second zone [x_VectorField_Z2, +inf]

	// Extrusion Parameters
	params_aero.nbr_couches = 4;			// Number of layer in extrusion
	params_aero.x_lim = 100;				// Limites physiques à partir desquelles
	params_aero.y_lim = -10000;			// l'insertion et la fusion de blocs
	params_aero.z_lim = -10000;			// sont autorisées

	// Smoothing Parameters
	params_aero.nbr_iter_smoothing_yao = 200 ;		// Number of iterations for the Yao Smoothing
	params_aero.damping_smoothing_yao = 0.2 ;		// Damping parameter for the Yao Smoothing

	// SU2 Writer Parameter
	params_aero.x_lim_SU2_inoutlet = -pow(10,6);		// Limit between inlet and outlet for SU2 writer

	// Mesh Generation
	AeroPipeline_2D algo_aero2D(params_aero);
	AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero2D_result);

}
*/

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
	//params_aero.input_file=dir+"/Aero/3D/biconique.vtk";
	params_aero.input_file=dir+"/Aero/3D/C2_3D_0.3.vtk";
	params_aero.output_file="AeroPipeline3D_Hexa.vtk";
	params_aero.output_dir="gmds/claire/tst/";

	// Extrusion Parameters
	params_aero.nbr_couches = 4;			// Number of layer in extrusion



	//---------------------//
	//    AERO PIPELINE    //
	//---------------------//
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