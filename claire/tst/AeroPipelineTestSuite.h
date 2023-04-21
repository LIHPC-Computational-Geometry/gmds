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
#include <gmds/blocking/CGNSWriter.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*                     CAS TEST 2D CLASSE AeroPipeline                        */
/*----------------------------------------------------------------------------*/

TEST(AeroPipelineTestClass, AeroPipeline2D_Test1)
{
	std::string dir(TEST_SAMPLES_DIR);
	std::string input_file=dir+"/Aero/2D/param_Apollo_2D.ini";

	// Mesh Generation
	AeroPipeline_2D algo_aero2D(input_file, dir);
	AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();


	blocking::CGNSWriter writer(algo_aero2D.getBlocking());
	writer.write("AeroPipeline_2D.cgns", "");

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero2D_result);

}

/*----------------------------------------------------------------------------*/
/*                    CAS TEST 3D CLASSE AeroPipeline                         */
/*----------------------------------------------------------------------------*/

TEST(AeroPipelineTestClass, AeroPipeline3D_Test1)
{
	std::string dir(TEST_SAMPLES_DIR);
	std::string input_file=dir+"/Aero/3D/param_C9_3D.ini";

	//---------------------//
	//    AERO PIPELINE    //
	//---------------------//
	AeroPipeline_3D algo_aero3D(input_file, dir);
	AbstractAeroPipeline::STATUS aero3D_result = algo_aero3D.execute();

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero3D_result);
}


/*----------------------------------------------------------------------------*/
/*            CAS TEST 3D CLASSE AeroPipeline avec Exceptions                 */
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