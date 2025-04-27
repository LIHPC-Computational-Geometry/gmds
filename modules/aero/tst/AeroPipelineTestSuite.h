//
// Created by rochec on 10/02/2022.
//

#include <gmds/aero/AbstractAeroPipeline.h>
#include <gmds/aero/AeroPipeline_2D.h>
#include <gmds/aero/AeroPipeline_3D.h>
#include <gmds/aero/AeroExtrusion_2D.h>
#include <gmds/aero/AeroException.h>
#include <gmds/aero/Params.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
#ifdef USE_CGNS
	#include <gmds/blocking/CGNSWriter.h>
#endif
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*                     CAS TEST 2D CLASSE AeroPipeline                        */
/*----------------------------------------------------------------------------*/

TEST(AeroPipelineTestClass, DISABLED_AeroPipeline2D_Test1)
{
	std::string dir(TEST_SAMPLES_DIR);
	//std::string input_file=dir+"/Aero/2D/param_Apollo_2D.ini";
	//std::string input_file=dir+"/Aero/2D/param_RAMCII_2D.ini";
	//std::string input_file=dir+"/Aero/2D/param_NACA_2D.ini";
	//std::string input_file=dir+"/Aero/2D/param_Orex_2D.ini";
	//std::string input_file=dir+"/Aero/2D/param_Stardust_2D_TEST.ini";
	//std::string input_file=dir+"/Aero/2D/param_Stardust_2D_param1.ini";
	//std::string input_file=dir+"/Aero/2D/param_Stardust_2D_v4.ini";
	std::string input_file=dir+"/Aero/2D/param_Diamond_Airfoil_2D.ini";
	//std::string input_file=dir+"/Aero/2D/param_U_Shaped.ini";

	// Mesh Generation
	AeroPipeline_2D algo_aero2D(input_file, dir);
	AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();

	#ifdef USE_CGNS
		blocking::CGNSWriter writer(algo_aero2D.getBlocking());
		writer.write("AeroPipeline_2D.cgns", "");
	#endif

	ASSERT_EQ(AbstractAeroPipeline::SUCCESS, aero2D_result);

}

/*----------------------------------------------------------------------------*/
/*                    CAS TEST 3D CLASSE AeroPipeline                         */
/*----------------------------------------------------------------------------*/

TEST(AeroPipelineTestClass, DISABLED_AeroPipeline3D_Test1)
{
	std::string dir(TEST_SAMPLES_DIR);
	//std::string input_file=dir+"/Aero/3D/param_Modified_CCF_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_Stardust_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_Apollo_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_Double_Ellipsoid_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_HyTRV_3D.ini";
	std::string input_file=dir+"/Aero/3D/param_RAMCII_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_Caretwing_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_C8_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_Ailerons_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_Ailerons3_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_TintinRocket_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_HiFIRE5_3D.ini";
	//std::string input_file=dir+"/Aero/3D/param_HB2_3D.ini";

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