/*---------------------------------------------------------------------------*/
#include <iostream>
#include <memory>
#include "gmds/rlBlocking/PolyCut.h"
#include "gmds/rlBlocking/MCTSAgentSPAM.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <unit_test_config.h>
/*---------------------------------------------------------------------------*/

std::string PATH_FOLDER = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_fix_cut_sheet/saveResults/results_PolyCut/";
/*---------------------------------------------------------------------------*/
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
int main() {
	gmds::cad::FACManager geom_model;
	std::string nameM= "cb2";
	set_up_SPAM(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);

	std::cout<<"list points geom size : "<<bl.geom_model()->getPoints().size()<<std::endl;

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	std::vector<double>  hist_empty;
	auto s = std::make_shared<PolyCutState>(&geom_model,&bl,hist_empty);

	PolyCutRewardFunction reward_function;
	MCTSAgent agent(&reward_function,1000000,1000);
	agent.run(s);

	std::cout<<"Nb runs: "<<agent.get_nb_iterations()-1,
	   std::cout<<", timing: "<<agent.get_nb_seconds()<<" s."<<std::endl;
	std::cout<<"==================== END TEST ! ===================="<<std::endl;
	auto best = std::dynamic_pointer_cast<PolyCutState> (agent.get_best_solution_uct());
	auto best_node = agent.get_best_node_uct();
	auto current_node =best_node;
	unsigned int numSave =0;
	while(current_node != nullptr){
		auto action =  std::dynamic_pointer_cast<PolyCutAction>(current_node->get_action());
		std::cout<<"Action: "<<std::endl;
		if(action!= nullptr) action->print();
		std::cout<<"Reward: "<<current_node->cumulative_reward<<std::endl;
		std::cout<<"Visits: "<<current_node->number_visits<<std::endl;
		std::cout<<"Nb children:"<<current_node->get_children().size()<<std::endl;
		auto current_state = std::dynamic_pointer_cast<PolyCutState> (current_node->get_state());
		current_state->m_blocking->save_vtk_blocking(std::to_string(numSave) + "OutPutCb0");
		numSave++;
		current_node = current_node->get_parent();
	}
	best->m_blocking->save_vtk_blocking("bestOutputCb0");
}