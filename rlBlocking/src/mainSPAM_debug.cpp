/*---------------------------------------------------------------------------*/
#include <iostream>
#include <memory>
#include "gmds/rlBlocking/PolyCut.h"
#include "gmds/rlBlocking/MCTSAgentSPAM.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <unit_test_config.h>
/*---------------------------------------------------------------------------*/

std::string PATH_FOLDER = "/home/legoffn/travail/gmds/build_20240321/output_paul";
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
	std::string nameM= "M1_ext";
	set_up_SPAM(&geom_model,nameM+".vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model, true);
//	std::string blockingInput = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_fix_cut_sheet/gmds/cmake-build-debug/bin/oldBis/inputCb2cut1.vtk";
//	gmds::blocking::CurvedBlocking bl(&geom_model, false);
//
//	//==================================================================
//	   // MESH READING
//	   //==================================================================
//	   std::cout<<"> Start mesh reading"<<std::endl;
//	//the used model is specified according to the class requirements.
//	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));
//
//	gmds::IGMeshIOService ioService2(&m);
//	gmds::VTKReader vtkReader2(&ioService2);
//	vtkReader2.setCellOptions(gmds::N|gmds::R);
//	vtkReader2.read(blockingInput);
//	gmds::MeshDoctor doc2(&m);
//	doc2.updateUpwardConnectivity();
//
//	//==================================================================
//
//	bl.init_from_mesh(m);

	std::cout<<"list points geom size : "<<bl.geom_model()->getPoints().size()<<std::endl;

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	std::cout<<"==================== BEGIN TEST : ===================="<<std::endl;

	std::vector<double>  hist_empty;
	auto s = std::make_shared<PolyCutState>(&geom_model,&bl,hist_empty);

	PolyCutRewardFunction reward_function;
	MCTSAgent agent(&reward_function,1000000,100000);
//	agent.activate_debug_mode("cb2_debug",MCTSAgent::OUT_ITERATION,1);
//	agent.run(s);

	auto m_tree  = new MCTSTree(s);
	s->m_blocking->save_vtk_blocking("s_init",0);

	auto actions = s->get_actions();
	std::cout<<"actions.size() "<<actions.size()<<std::endl;
	for(auto a: actions) {
//		std::cout<< std::dynamic_pointer_cast<PolyCutAction>(a)->print()<<std::endl;
		std::dynamic_pointer_cast<PolyCutAction>(a)->print();
	}

	auto children = m_tree->get_children();
	std::cout<<"children.size() "<<children.size()<<std::endl;

	std::cout<<"m_tree->is_fully_expanded() "<<m_tree->is_fully_expanded()<<std::endl;
	std::cout<<"m_tree->is_terminal() "<<m_tree->is_terminal()<<std::endl;

	m_tree->expand();
	children = m_tree->get_children();
	std::cout<<"children.size() "<<children.size()<<std::endl;

	auto child = children[0];
	auto s_c = child->get_state();
	std::dynamic_pointer_cast<PolyCutState>(s_c)->m_blocking->save_vtk_blocking("s_c",0);

	auto actions_c = s_c->get_actions();
	for(auto a: actions_c) {
		std::dynamic_pointer_cast<PolyCutAction>(a)->print();
	}

	child->expand();
	children = child->get_children();
	std::cout<<"children.size() "<<children.size()<<std::endl;
	child = children[0];
	s_c = child->get_state();
	std::dynamic_pointer_cast<PolyCutState>(s_c)->m_blocking->save_vtk_blocking("s_cc",0);

//	auto node = agent.select(m_tree);

//	auto node = select(m_tree);
//	// 2. EXPAND by adding a single child (if not terminal or not fully expanded)
//	node = expand(node);
//	// 3. SIMULATE (if not terminal)
//	auto reward = simulate(node);
//	// 4. BACK PROPAGATION
//	back_propagate(node,reward);

//	std::cout<<"Nb runs: "<<agent.get_nb_iterations()-1,
//	   std::cout<<", timing: "<<agent.get_nb_seconds()<<" s."<<std::endl;
//	std::cout<<"==================== END TEST ! ===================="<<std::endl;
//	auto best = std::dynamic_pointer_cast<PolyCutState> (agent.get_best_solution_uct());
//	auto best_node = agent.get_best_node_uct();
//	auto current_node =best_node;
//	unsigned int numSave =0;
//	while(current_node != nullptr){
//		auto action =  std::dynamic_pointer_cast<PolyCutAction>(current_node->get_action());
//		std::cout<<"Action: "<<std::endl;
//		if(action!= nullptr) action->print();
//		std::cout<<"Reward: "<<current_node->cumulative_reward<<std::endl;
//		std::cout<<"Visits: "<<current_node->number_visits<<std::endl;
//		std::cout<<"Nb children:"<<current_node->get_children().size()<<std::endl;
//		auto current_state = std::dynamic_pointer_cast<PolyCutState> (current_node->get_state());
//		//current_state->m_blocking->save_vtk_blocking(std::to_string(numSave) + "OutPutCb0");
//		numSave++;
//		current_node = current_node->get_parent();
//	}
//	best->m_blocking->save_vtk_blocking("bestOutputCb0");
}