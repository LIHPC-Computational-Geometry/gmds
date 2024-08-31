/*---------------------------------------------------------------------------*/
#include <mcts/MCTSLoop.h>
#include <iostream>
/*---------------------------------------------------------------------------*/
MCTSLoop::MCTSLoop(MCTSAgent &AAgent,
                   std::shared_ptr<IState> ARootState,
                   MCTSLoop::NODE_SELECTION ASelectionMode,
                   const int AMaxIterations,
                   const bool ADisplayInfo)
: m_agent(AAgent),
  m_root_state(ARootState),
  m_selection_mode(ASelectionMode),
  m_max_iterations(AMaxIterations),
  m_display_info(ADisplayInfo),
  m_nb_iterations(0)
{}
/*---------------------------------------------------------------------------*/
MCTSLoop::~MCTSLoop() {}
/*---------------------------------------------------------------------------*/
void MCTSLoop::run() {
    auto current_state = m_root_state;
    m_nb_iterations =0;
    for(auto i=0;i<m_max_iterations&& !current_state->win();i++) {
        m_agent.run(current_state);
        m_nb_iterations += m_agent.get_nb_iterations() - 1;
        if (m_display_info) {
            std::cout << "Iteration " << i << ", nb runs: " << m_agent.get_nb_iterations() - 1,
                    std::cout << ", timing: " << m_agent.get_nb_seconds() << " s." << std::endl;
        }
        if (m_selection_mode == MCTSLoop::BEST_SOLUTION)
            current_state = m_agent.get_best_solution();
        else {
			   current_state = m_agent.get_most_visited_child();
		  }
        if (m_display_info) {
            if (current_state->win()) {
                std::cout << "\t found a winning solution" << std::endl;
            }
        }
    }

}
/*---------------------------------------------------------------------------*/
