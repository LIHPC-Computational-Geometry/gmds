/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/BlockingClassifier.h>
#include <gmds/mctsblock/BlockingRewardFunction.h>
#include <gmds/mctsblock/BlockingState.h>
#include <mcts/MCTSAgent.h>
#include <mcts/MCTSSelectionFunction.h>
#include <gmds/mctsblock/Blocking.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
#include <nlohmann/json.hpp>
using json = nlohmann::json;
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
json read_param_file(const std::string& AFileName) {
	json data;
	// read in the json file
	std::ifstream json_file(AFileName);
	data = json::parse(json_file); // initialize json object with what was read from file

	json_file.close();

	return data;
}
/*----------------------------------------------------------------------------*/
void
init_geom(cad::FACManager &AGeomManager, const std::string& AGeomFile)
{
	Mesh m_vol(MeshModel(DIM3 | R | F | E | N | R2N | R2F | R2E | F2N | F2R | F2E
	                                 | E2F | E2N | N2E));
	IGMeshIOService ioService(&m_vol);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N | R);
	vtkReader.read(AGeomFile);
	MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	AGeomManager.initFrom3DMesh(&m_vol);
}
/*----------------------------------------------------------------------------*/
void
display_info(std::shared_ptr<mctsblock::BlockingState> state)
{
	auto bl = state->get_blocking();
	std::cout << "Blocking " << bl.get() << " (N,E,F,N): " << bl->get_nb_cells<0>() << ", " << bl->get_nb_cells<1>() << ", " << bl->get_nb_cells<2>() << ", "
	          << bl->get_nb_cells<3>() << std::endl;
	auto errors = mctsblock::BlockingClassifier(bl.get()).detect_classification_errors();
	std::cout << "\t non captured points (" << errors.non_captured_points.size() << "): ";
	for (auto i : errors.non_captured_points) {
		auto pi = bl->geom_model()->getPoint(i)->point();
		std::cout << i << " ("<<pi.X()<<", "<<pi.Y()<<","<<pi.Z()<<") ";
	}
	std::cout << std::endl;
	std::cout << "\t non captured curves (" << errors.non_captured_curves.size() << "): ";
	for (auto i : errors.non_captured_curves)
		std::cout << i << " ";
	std::cout << std::endl;
	std::cout << "\t non captured surfaces (" << errors.non_captured_surfaces.size() << "): ";
	for (auto i : errors.non_captured_surfaces)
		std::cout << i << " ";
	std::cout << std::endl;
	std::cout << "\t score: " << state->computeScore() << std::endl;
}
/*----------------------------------------------------------------------------*/
void
read_blocks(Mesh &AMesh, const std::string AFileName)
{
	IGMeshIOService ioService(&AMesh);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N | R);
	vtkReader.read(AFileName);
	MeshDoctor doc(&AMesh);
	doc.updateUpwardConnectivity();
}
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "============== MCTS Blocker ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
    std::string file_geom, file_block_out, file_block_in, file_param;
    if (argc != 4 && argc !=5) {
        std::cout << "Wrong number of parameters. It require 4 or 5 parameters : \n";
		  std::cout << "  1. [IN]\t MCTS parameter file (.json), \n";
		  std::cout << "  2. [IN]\t tetrahedral mesh (.vtk) that describes the geometry, \n";
		  std::cout << "  3. [OUT]\t the final blocks (.vtk),\n";
		  std::cout << "  4. [IN/OPT]\t an init block structure (.vtk)"<< std::endl;
		  exit(0);
    }
	 file_param = std::string(argv[1]);
	 file_geom = std::string(argv[2]);
    file_block_out = std::string(argv[3]);

	 auto with_init_blocking = false;
	 if(argc==5){
		  with_init_blocking = true;
		  file_block_in = std::string(argv[4]);
	 }
	 std::cout << "=======================================" << std::endl;
    std::cout << "Parameters " << std::endl;
	 std::cout << "  - Parameter file:\t " << file_param << std::endl;
	 std::cout << "  - Geometry file:\t " << file_geom << std::endl;
	 if(with_init_blocking)
		  std::cout << "  - Input block file:\t " << file_block_out << std::endl;
    std::cout << "  - Output block file:\t " << file_block_out << std::endl;
    std::cout << "=======================================" << std::endl;

	 //==================================================================
	 // INITIALIZATION

	 auto params = read_param_file(file_param);
	 std::cout<<"Nb MCTS iterations:\t\t "<<params.at("nb_mcts_iter")<<std::endl;
	 std::cout<<"Nb loop iterations:\t\t "<<params.at("nb_loop_iter")<<std::endl;
	 std::cout<<"UCT C value:\t\t\t "<<params.at("uct_C")<<std::endl;
	 std::cout<<"MCTS Simulation depth:\t\t "<<params.at("simulation_depth")<<std::endl;
	 std::cout<<"MCTS iteration max time (s):\t "<<params.at("max_mcts_iteration_time")<<std::endl;
	 std::cout << "=======================================" << std::endl;
	 
	 // Geometry initialization
	 cad::FACManager geom_model;
	 init_geom(geom_model, file_geom);
	 geom_model.buildGTSTree();
	 
	 // We initialize the blocking structure from the geom model
	 mctsblock::Blocking bl(&geom_model, false);
	 // and we add the initial block structure afterward
	 if(with_init_blocking) {
		  Mesh init_blocks(MeshModel(DIM3 | R | N | R2N | N2R));
		  read_blocks(init_blocks, file_block_in);
		  bl.init_from_mesh(init_blocks);
	 }
	 else{
		  //we use the bounding box
		  bl.init_from_bounding_box();
	 }


	 mctsblock::BlockingRewardFunction reward_function;

	 auto init_state = std::make_shared<mctsblock::BlockingState>(std::make_shared<mctsblock::Blocking>(bl));
	 auto optimal_score = init_state->get_expected_optimal_score();
	 std::cout << "Expected optimal score: " << optimal_score << std::endl;
	 std::cout << "=======================================" << std::endl;
	 SPUCTSelectionFunction select_function(params.at("uct_C"), optimal_score);

	 MCTSAgent agent(&reward_function, &select_function,
	                 params.at("nb_mcts_iter"),
	                 params.at("max_mcts_iteration_time"),
	                 params.at("simulation_depth"));

	 agent.activate_debug_mode("blocking", MCTSAgent::OUT_END_ONLY, 1000);
	 auto current_state = init_state;
	 display_info(init_state);
	 auto nb_loop_iter = params.at("nb_loop_iter");

	 for (auto i = 0; i < nb_loop_iter && !current_state->win() && !current_state->lost(); i++) {
		  agent.run(current_state);
		  std::cout << "=======================================" << std::endl;

		  std::cout << "Iteration " << i << ", nb runs: " << agent.get_nb_iterations() - 1;
		  std::cout << ", timing: " << agent.get_nb_seconds() << " s." << std::endl;

		  current_state = std::dynamic_pointer_cast<mctsblock::BlockingState>(agent.get_most_winning_child());
		  display_info(current_state);
		  current_state->get_blocking()->save_vtk_blocking("loop_" + std::to_string(i));

	 }
	 std::cout << "=======================================" << std::endl;

	 if (init_state->win()) {
		  std::cout << "Init is a best solution : " << init_state->computeScore() << std::endl;
		  init_state->get_blocking()->save_vtk_blocking(file_block_out);
	 }
	 else {
		  auto solution = *std::dynamic_pointer_cast<mctsblock::BlockingState>(agent.get_best_winning_solution());
		  std::cout << "Best solution score: " << solution.computeScore() <<" (optimal: "<<optimal_score<<")"<<std::endl;
		  solution.get_blocking()->save_vtk_blocking(file_block_out);
	 }
	 if(current_state->win()) {
		 std::cout << "\n\t >>>>>> WIN <<<<<<" << std::endl;
		 return 1;
	 }
	 else if(current_state->lost()) {
		 std::cout << "\n\t >>>>>> LOST <<<<<<" << std::endl;
		 return 0;
	 }

}