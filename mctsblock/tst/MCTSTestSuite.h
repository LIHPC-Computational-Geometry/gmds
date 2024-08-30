/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/mctsblock/Blocking.h>
#include <gmds/mctsblock/BlockingClassifier.h>
#include <unit_test_config.h>

#include <mcts/MCTSAgent.h>
#include <mcts/MCTSLoop.h>
#include <mcts/MCTSSelectionFunction.h>
#include <mcts/MCTSTree.h>

#include <gmds/mctsblock/BlockingAction.h>
#include <gmds/mctsblock/BlockingRewardFunction.h>
#include <gmds/mctsblock/BlockingState.h>
/*----------------------------------------------------------------------------*/
void
setUpMCTS(gmds::cad::FACManager &AGeomManager)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
	                                 | gmds::E2F | gmds::E2N | gmds::N2E));
	//	std::string dir(TEST_SAMPLES_DIR);
	//
	//
	//	std::string vtk_file = dir + "/tet_in_box.vtk";
	//	std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/cube.vtk";
	//	std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/cube_minus_edge.vtk";
	// std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/cube_minus_corner.vtk";
	// std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/cube_minus_two_edges.vtk";
	// std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/U_shape.vtk";
	// std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/XYZ.vtk";
	// std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/T_shape.vtk";

	//	std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/hole.vtk";
	std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/XXYZ.vtk";
	//	std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/cross_3D.vtk";

	gmds::IGMeshIOService ioService(&m_vol);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	gmds::MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	AGeomManager.initFrom3DMesh(&m_vol);
}
/*----------------------------------------------------------------------------*/
std::tuple<int, int, int, int>
stat(gmds::mctsblock::Blocking &ABlocking)
{
	auto nb_on_vertex = 0;
	auto nb_on_curve = 0;
	auto nb_on_surface = 0;
	auto nb_in_volume = 0;
	std::vector<gmds::mctsblock::Blocking::Node> all_nodes = ABlocking.get_all_nodes();
	for (auto n : all_nodes) {
		if (n->info().geom_dim == 0)
			nb_on_vertex++;
		else if (n->info().geom_dim == 1)
			nb_on_curve++;
		else if (n->info().geom_dim == 2)
			nb_on_surface++;
		else if (n->info().geom_dim == 3)
			nb_in_volume++;
	}
	return std::make_tuple(nb_on_vertex, nb_on_curve, nb_on_surface, nb_in_volume);
}
/*----------------------------------------------------------------------------*/
void
display_info(std::shared_ptr<gmds::mctsblock::BlockingState> state)
{
	auto bl = state->get_blocking();
	std::cout << "Blocking " << bl.get() << " (N,E,F,N): " << bl->get_nb_cells<0>() << ", " << bl->get_nb_cells<1>() << ", " << bl->get_nb_cells<2>() << ", "
	          << bl->get_nb_cells<3>() << std::endl;
	auto errors = gmds::mctsblock::BlockingClassifier(bl.get()).detect_classification_errors();
	std::cout << "\t non captured points (" << errors.non_captured_points.size() << "): ";
	for (auto i : errors.non_captured_points)
		std::cout << i << " ";
	std::cout << std::endl;
	std::cout << "\t non captured curves (" << errors.non_captured_curves.size() << "): ";
	for (auto i : errors.non_captured_curves)
		std::cout << i << " ";
	std::cout << std::endl;
	std::cout << "\t score: " << state->computeScore() << std::endl;
}
/*----------------------------------------------------------------------------*/

void
read_blocks(gmds::Mesh &AMesh, const std::string AFileName)
{
	gmds::IGMeshIOService ioService(&AMesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(AFileName);
	gmds::MeshDoctor doc(&AMesh);
	doc.updateUpwardConnectivity();
}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, box)
{
	gmds::cad::FACManager geom_model;
	setUpMCTS(geom_model);
	gmds::mctsblock::Blocking bl(&geom_model, true);

	/*		gmds::Mesh init_blocks(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::R2N | gmds::N2R));
	      read_blocks(init_blocks, "/Users/bibi/Developpement/MCTS-data/solutions/solution_hole.vtk");
	      bl.init_from_mesh(init_blocks);
	   */
	ASSERT_EQ(8, bl.get_all_nodes().size());
	ASSERT_EQ(12, bl.get_all_edges().size());
	ASSERT_EQ(6, bl.get_all_faces().size());
	ASSERT_EQ(1, bl.get_all_blocks().size());
	auto nb_geom_points = geom_model.getNbPoints();
	auto nb_geom_curves = geom_model.getNbCurves();

	auto optimal_score = gmds::mctsblock::BlockingState::weight_nodes * nb_geom_points + gmds::mctsblock::BlockingState::weight_edges * nb_geom_curves;
	std::cout << "Expected optimal score: " << optimal_score << std::endl;
	auto nb_mcts_iter = 1000;
	auto nb_loop_iter = 100;
	bool get_best_solution = false;
	gmds::mctsblock::BlockingRewardFunction reward_function;
	//UCBSelectionFunction select_function;
	SPUCTSelectionFunction select_function(1.42, optimal_score);
	auto init_state = std::make_shared<gmds::mctsblock::BlockingState>(std::make_shared<gmds::mctsblock::Blocking>(bl));

	MCTSAgent agent(&reward_function, &select_function, nb_mcts_iter, 100, 100);

	agent.activate_debug_mode("blocking", MCTSAgent::OUT_END_ONLY, 1000);
	auto current_state = init_state;
	display_info(init_state);
	for (auto i = 0; i < nb_loop_iter && !current_state->win() && !current_state->lost(); i++) {
		agent.run(current_state);

		std::cout << "Iteration " << i << ", nb runs: " << agent.get_nb_iterations() - 1;
		std::cout << ", timing: " << agent.get_nb_seconds() << " s." << std::endl;

		current_state = std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_most_visited_child());
		display_info(current_state);
		current_state->get_blocking()->save_vtk_blocking("loop_" + std::to_string(i));
		if (current_state->win()) {
			std::cout << "\t found a winning solution" << std::endl;
			auto mem = current_state->get_memory();
			std::cout << "Memory: " << std::endl;
			for (auto m : mem)
				std::cout << m << " ";
			std::cout << std::endl;
			std::cout << "Score: " << current_state->computeScore() << std::endl;
		}
	}
	if (init_state->win()) {
		std::cout << "Init is a best solution : " << init_state->computeScore() << std::endl;
		init_state->get_blocking()->save_vtk_blocking("best_init.vtk");
	}
	else {
		auto solution = *std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_best_solution());
		std::cout << "Best solution score: " << solution.computeScore() << std::endl;
		solution.get_blocking()->save_vtk_blocking("best.vtk");
	}
	if(current_state->win())
		std::cout<<"Final WIN"<<std::endl;
	else if(current_state->lost())
		std::cout<<"Final LOST"<<std::endl;
}