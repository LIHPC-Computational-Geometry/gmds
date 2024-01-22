#ifndef GMDS_MCTSTESTSUITE_H
#define GMDS_MCTSTESTSUITE_H
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSAlgorithm.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
/**@brief setup function that initialize a geometric model using the faceted
 * representation and an input vtk file name. The vtk file must contain a
 * tetrahedral mesh
 *
 * @param AGeomModel geometric model we initialize
 * @param AFileName vtk filename
 */
void set_up_MCTS(gmds::cad::FACManager* AGeomModel, const std::string AFileName)
{
	gmds::Mesh vol_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
	                                    | gmds::E2F | gmds::E2N | gmds::N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir +"/"+ AFileName;
	gmds::IGMeshIOService ioService(&vol_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	gmds::MeshDoctor doc(&vol_mesh);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	AGeomModel->initFrom3DMesh(&vol_mesh);

}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, testExAgloCb2)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "cb2";
	set_up_MCTS(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking("/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/cb2/cb2_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	MCTSAlgorithm *algo = new MCTSAlgorithm(&geom_model,&bl);
	algo->execute(nameM);


	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}

/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, testExAgloM1)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "M1";
	set_up_MCTS(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking("/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/M1/M1_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	MCTSAlgorithm *algo = new MCTSAlgorithm(&geom_model,&bl);
	algo->execute(nameM);


	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, testExAgloCb1)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "cb1";
	set_up_MCTS(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking("/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/cb1/cb1_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	MCTSAlgorithm *algo = new MCTSAlgorithm(&geom_model,&bl);
	algo->execute(nameM);


	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, testExAgloCb3)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "cb3";
	set_up_MCTS(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking("/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/cb3/cb3_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	MCTSAlgorithm *algo = new MCTSAlgorithm(&geom_model,&bl);
	algo->execute(nameM);


	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, testExAgloCb4)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "cb4";
	set_up_MCTS(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking("/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/cb4/cb4_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	MCTSAlgorithm *algo = new MCTSAlgorithm(&geom_model,&bl);
	algo->execute(nameM);


	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, testExAgloCb5)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "cb5";
	set_up_MCTS(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking("/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/cb5/cb5_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	MCTSAlgorithm *algo = new MCTSAlgorithm(&geom_model,&bl);
	algo->execute(nameM);


	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}

/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, testExAgloB0)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "B0";
	set_up_MCTS(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking("/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/B0/B0_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	MCTSAlgorithm *algo = new MCTSAlgorithm(&geom_model,&bl);
	algo->execute(nameM);


	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_MCTSTESTSUITE_H
