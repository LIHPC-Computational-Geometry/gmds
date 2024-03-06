#ifndef GMDS_SPAMTESTSUITE_H
#define GMDS_SPAMTESTSUITE_H
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "gmds/rlBlocking/MCTSAgentSPAM.h"
#include "gmds/rlBlocking/PolyCut.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <unit_test_config.h>

std::string PATH_FOLDER = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_fix_cut_sheet/saveResults/results_PolyCut/";
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
void set_up_SPAM(gmds::cad::FACManager* AGeomModel, const std::string AFileName)
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
TEST(SPAMTestSuite, testSPAMcb2)
{

	gmds::cad::FACManager geom_model;
	std::string nameM= "cb2";
	set_up_SPAM(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	bl.save_vtk_blocking(PATH_FOLDER+"cb2/cb2_init_blocking.vtk");

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	std::vector<double>  hist_empty;
	auto s = std::make_shared<PolyCutState>(&geom_model,&bl,hist_empty);

	PolyCutRewardFunction reward_function;
	MCTSAgent agent(&reward_function,10000);
	agent.run(s);

	std::cout<<"Nb runs: "<<agent.get_nb_iterations()-1,
	   std::cout<<", timing: "<<agent.get_nb_seconds()<<" s."<<std::endl;
	std::cout<<"==================== END TEST ! ===================="<<std::endl;


}
#endif     // GMDS_SPAMTESTSUITE_H
