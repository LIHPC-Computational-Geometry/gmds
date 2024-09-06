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
#include <mctsblock/config.h>
/*----------------------------------------------------------------------------*/
void
init_geom(gmds::cad::FACManager &AGeomManager, const std::string& AFileName)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
	                                 | gmds::E2F | gmds::E2N | gmds::N2E));
	std::string dir(MCTS_DATA_DIR);

	gmds::IGMeshIOService ioService(&m_vol);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(dir + AFileName);
	gmds::MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	AGeomManager.initFrom3DMesh(&m_vol);
	AGeomManager.write_surfaces("geometry");
}
/*----------------------------------------------------------------------------*/
void read_blocks(gmds::Mesh &AMesh, const std::string AFileName)
{
	gmds::IGMeshIOService ioService(&AMesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(AFileName);
	gmds::MeshDoctor doc(&AMesh);
	doc.updateUpwardConnectivity();
}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, cube) {
	gmds::cad::FACManager geom_model;
	init_geom(geom_model, "cube.vtk");
	gmds::mctsblock::Blocking bl(&geom_model, true);
	auto nb_mcts_iter = 1000;
	auto nb_loop_iter = 100;


	auto init_state = std::make_shared<gmds::mctsblock::BlockingState>(std::make_shared<gmds::mctsblock::Blocking>(bl));
	auto optimal_score = init_state->get_expected_optimal_score();
	std::cout<<"Cube exp. score = "<<optimal_score<<std::endl;

	SPUCTSelectionFunction select_function(1.42, optimal_score);
	gmds::mctsblock::BlockingRewardFunction reward_function;

	MCTSAgent agent(&reward_function, &select_function, nb_mcts_iter, 100, 20);
	auto current_state = init_state;
	for (auto i = 0; i < nb_loop_iter && !current_state->win() && !current_state->lost(); i++) {
		agent.run(current_state);
		current_state = std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_most_winning_child());
	}
	std::shared_ptr<gmds::mctsblock::Blocking> final_blocking;
	if (init_state->win()) {
		final_blocking = init_state->get_blocking();
	}
	else {
		auto solution = *std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_best_winning_solution());
		final_blocking = solution.get_blocking();
	}

	ASSERT_EQ(final_blocking->get_nb_cells<3>(),1);
	ASSERT_TRUE(current_state->win());
}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, cube_minus_edge) {
	gmds::cad::FACManager geom_model;
	init_geom(geom_model, "cube_minus_edge.vtk");
	gmds::mctsblock::Blocking bl(&geom_model, true);

	auto nb_mcts_iter = 1000;
	auto nb_loop_iter = 100;
	gmds::mctsblock::BlockingRewardFunction reward_function;

	auto init_state = std::make_shared<gmds::mctsblock::BlockingState>(std::make_shared<gmds::mctsblock::Blocking>(bl));
	auto optimal_score = init_state->get_expected_optimal_score();
	SPUCTSelectionFunction select_function(1.42, optimal_score);

	MCTSAgent agent(&reward_function, &select_function, nb_mcts_iter, 100, 20);
	auto current_state = init_state;
	for (auto i = 0; i < nb_loop_iter && !current_state->win() && !current_state->lost(); i++) {
		agent.run(current_state);
		current_state = std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_most_winning_child());
	}
	std::shared_ptr<gmds::mctsblock::Blocking> final_blocking;
	if (init_state->win()) {
		final_blocking = init_state->get_blocking();
	}
	else {
		auto solution = *std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_best_winning_solution());
		final_blocking = solution.get_blocking();
	}

	ASSERT_EQ(final_blocking->get_nb_cells<3>(),3);
	ASSERT_TRUE(current_state->win());
}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, cube_minus_corner) {
	gmds::cad::FACManager geom_model;
	init_geom(geom_model, "cube_minus_corner.vtk");
	gmds::mctsblock::Blocking bl(&geom_model, true);

	auto nb_mcts_iter = 1000;
	auto nb_loop_iter = 100;
	gmds::mctsblock::BlockingRewardFunction reward_function;
	auto init_state = std::make_shared<gmds::mctsblock::BlockingState>(std::make_shared<gmds::mctsblock::Blocking>(bl));
	auto optimal_score = init_state->get_expected_optimal_score();
	SPUCTSelectionFunction select_function(1.42, optimal_score);

	MCTSAgent agent(&reward_function, &select_function, nb_mcts_iter, 100, 20);
	auto current_state = init_state;
	for (auto i = 0; i < nb_loop_iter && !current_state->win() && !current_state->lost(); i++) {
		agent.run(current_state);
		current_state = std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_most_winning_child());
	}
	std::shared_ptr<gmds::mctsblock::Blocking> final_blocking;
	if (init_state->win()) {
		final_blocking = init_state->get_blocking();
	}
	else {
		auto solution = *std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_best_winning_solution());
		final_blocking = solution.get_blocking();
	}

	ASSERT_EQ(final_blocking->get_nb_cells<3>(),7);
	ASSERT_TRUE(current_state->win());
}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, cube_minus_two_edges) {
	gmds::cad::FACManager geom_model;
	init_geom(geom_model, "cube_minus_two_edges.vtk");
	gmds::mctsblock::Blocking bl(&geom_model, true);

	auto nb_mcts_iter = 1000;
	auto nb_loop_iter = 100;
	gmds::mctsblock::BlockingRewardFunction reward_function;
	auto init_state = std::make_shared<gmds::mctsblock::BlockingState>(std::make_shared<gmds::mctsblock::Blocking>(bl));
	auto optimal_score = init_state->get_expected_optimal_score();
	SPUCTSelectionFunction select_function(1.42, optimal_score);

	MCTSAgent agent(&reward_function, &select_function, nb_mcts_iter, 100, 20);
	auto current_state = init_state;
	for (auto i = 0; i < nb_loop_iter && !current_state->win() && !current_state->lost(); i++) {
		agent.run(current_state);
		current_state = std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_most_winning_child());
	}
	std::shared_ptr<gmds::mctsblock::Blocking> final_blocking;
	if (init_state->win()) {
		final_blocking = init_state->get_blocking();
	}
	else {
		auto solution = *std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_best_winning_solution());
		final_blocking = solution.get_blocking();
	}
	ASSERT_EQ(final_blocking->get_nb_cells<3>(),5);
	ASSERT_TRUE(current_state->win());
}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, U_shape) {
	gmds::cad::FACManager geom_model;
	init_geom(geom_model, "U_shape.vtk");
	gmds::mctsblock::Blocking bl(&geom_model, true);

	auto nb_mcts_iter = 1000;
	auto nb_loop_iter = 100;
	gmds::mctsblock::BlockingRewardFunction reward_function;
	auto init_state = std::make_shared<gmds::mctsblock::BlockingState>(std::make_shared<gmds::mctsblock::Blocking>(bl));
	auto optimal_score = init_state->get_expected_optimal_score();
	SPUCTSelectionFunction select_function(1.42, optimal_score);

	MCTSAgent agent(&reward_function, &select_function, nb_mcts_iter, 100, 20);
	auto current_state = init_state;
	for (auto i = 0; i < nb_loop_iter && !current_state->win() && !current_state->lost(); i++) {
		agent.run(current_state);
		current_state = std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_most_winning_child());
	}
	std::shared_ptr<gmds::mctsblock::Blocking> final_blocking;
	if (init_state->win()) {
		final_blocking = init_state->get_blocking();
	}
	else {
		auto solution = *std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_best_winning_solution());
		final_blocking = solution.get_blocking();
	}
	ASSERT_EQ(final_blocking->get_nb_cells<3>(),5);
	ASSERT_TRUE(current_state->win());
}
/*----------------------------------------------------------------------------*/
//std::string vtk_file = dir + "XYZ.vtk";
//std::string vtk_file = dir + "T_shape.vtk";
//std::string vtk_file = dir + "hole.vtk";
//std::string vtk_file = dir + "XXYZ.vtk";
//std::string vtk_file = dir + "cross_3D.vtk";
//std::string vtk_file = dir + "cylinder_quart.vtk";
//std::string vtk_file = dir + "curved_shape_1.vtk"
//std::string vtk_file = dir + "cube_minus_cyl_edge.vtk";
//std::string vtk_file = dir + "B0.vtk";