/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/Blocking.h>
#include <gmds/mctsblock/BlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>

#include <mcts/MCTSAgent.h>
#include <mcts/MCTSTree.h>
#include <mcts/MCTSLoop.h>
#include <mcts/MCTSSelectionFunction.h>

#include <gmds/mctsblock/BlockingState.h>
#include <gmds/mctsblock/BlockingRewardFunction.h>
#include <gmds/mctsblock/BlockingAction.h>
/*----------------------------------------------------------------------------*/
void
setUpMCTS(gmds::cad::FACManager &AGeomManager)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N |
	                                 gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N |
	                                 gmds::F2R | gmds::F2E
	                                 | gmds::E2F | gmds::E2N | gmds::N2E));
//	std::string dir(TEST_SAMPLES_DIR);
//	std::string vtk_file = dir + "/tet_in_box.vtk";
	std::string vtk_file = "/Users/bibi/Developpement/MCTS-data/C1.vtk";

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
std::tuple<int,int,int,int> stat(gmds::mctsblock::Blocking& ABlocking){
    auto nb_on_vertex=0;
    auto nb_on_curve=0;
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gmds::mctsblock::Blocking::Node> all_nodes = ABlocking.get_all_nodes();
    for(auto n:all_nodes){
        if(n->info().geom_dim==0)
            nb_on_vertex++;
        else if(n->info().geom_dim==1)
            nb_on_curve++;
        else if(n->info().geom_dim==2)
            nb_on_surface++;
        else if(n->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_vertex,nb_on_curve,nb_on_surface,nb_in_volume);
}
/*----------------------------------------------------------------------------*/
TEST(MCTSTestSuite, box)
{
	 gmds::cad::FACManager geom_model;
	 setUpMCTS(geom_model);
	 gmds::mctsblock::Blocking bl(&geom_model, true);

	 ASSERT_EQ(8, bl.get_all_nodes().size());
	 ASSERT_EQ(12,bl.get_all_edges().size());
	 ASSERT_EQ(6, bl.get_all_faces().size());
	 ASSERT_EQ(1, bl.get_all_blocks().size());


	 auto nb_mcts_iter = 1000;
	 auto nb_loop_iter = 100;
	 bool get_best_solution = false;
	 gmds::mctsblock::BlockingRewardFunction reward_function;
	 UCBSelectionFunction select_function;
	 auto init_state = std::make_shared<gmds::mctsblock::BlockingState>(&bl);

	 MCTSAgent agent(&reward_function, &select_function, nb_mcts_iter);
	 agent.activate_debug_mode("blocking", MCTSAgent::OUT_END_ONLY, 1000);

	 MCTSLoop loop(agent, init_state, MCTSLoop::BEST_CHILD, nb_loop_iter, true);
	 loop.run();

	 auto solution =  *std::dynamic_pointer_cast<gmds::mctsblock::BlockingState>(agent.get_best_solution());
	 std::cout<<"Best solution score: "<<solution.computeScore()<<std::endl;
	 solution.get_blocking()->save_vtk_blocking("toto.vtk");

}