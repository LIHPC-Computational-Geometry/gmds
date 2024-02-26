/*---------------------------------------------------------------------------*/
#include "gmds/rlBlocking/MCTSAgentSPAM.h"
#include <chrono>
#include <iostream>
#include <random>
/*---------------------------------------------------------------------------*/
MCTSAgent::MCTSAgent(const IRewardFunction* ARewardFunction,
                     const int AMaxIterations,
                     const int AMaxSeconds,
                     const int AMaxSimulationDepth)
: m_tree(nullptr),
m_reward_function(ARewardFunction),
m_max_iterations(AMaxIterations),
m_max_seconds(AMaxSeconds),
m_simulation_depth(AMaxSimulationDepth)
{}
/*---------------------------------------------------------------------------*/
MCTSAgent::~MCTSAgent() {
    delete m_tree;
}
/*---------------------------------------------------------------------------*/
void MCTSAgent::run(std::shared_ptr<IState> ARootState) {
    //Build the initial tree
   m_tree  = new MCTSTree(ARootState);

    int i=0;
    auto time0 = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::ratio<1>> elapsed= std::chrono::steady_clock::now()-time0;

    while (i<=m_max_iterations && elapsed.count() <= m_max_seconds){
        // 1 SELECTION - we explore/exploit the existing tree to find a node to work on
        auto node = select(m_tree);
        // 2. EXPAND by adding a single child (if not terminal or not fully expanded)
        node = expand(node);
        // 3. SIMULATE (if not terminal)
        auto reward = simulate(node);
        // 4. BACK PROPAGATION
        back_propagate(node,reward);
        //increase loop counters
        i++;
        elapsed= std::chrono::steady_clock::now()-time0;
    }
    m_nb_iterations=i;
    m_nb_seconds =elapsed.count();
}
/*---------------------------------------------------------------------------*/
std::shared_ptr<IState> MCTSAgent::get_best_solution() {
    const MCTSTree* node = m_tree;
    while (!node->is_terminal() && node->has_children()){
        node=node->get_most_visited_child();
    }
    return node->get_state();
}
/*---------------------------------------------------------------------------*/
MCTSTree* MCTSAgent::select(MCTSTree* ANode) {
    MCTSTree *node = ANode;
    while (node->is_fully_expanded() && !node->is_terminal()){
        node = node->get_best_uct_child();
    }
    //We have reached a terminal node here
    return  node;

}
/*---------------------------------------------------------------------------*/
MCTSTree* MCTSAgent::expand(MCTSTree* ANode) {
    if(!ANode->is_fully_expanded() && !ANode->is_terminal())
        return  ANode->expand();

    if(ANode->is_terminal())
        return ANode;

    return nullptr;
}
/*---------------------------------------------------------------------------*/
double MCTSAgent::simulate(MCTSTree* ANode) {
    auto state =ANode->get_state();

    if(!ANode->is_terminal()) {
        for (int d = 0; d < m_simulation_depth; d++) {
            if (!state->is_terminal()) {
                auto a = get_random_action(state);
                state = state->apply(a);
            }
        }
    }
    return m_reward_function->evaluate(state);
}
/*---------------------------------------------------------------------------*/
void MCTSAgent::back_propagate(MCTSTree* ANode, double AReward) {
    auto  node = ANode;
    while(node){
        node->cumulative_reward+=AReward;
        node->number_visits+=1;
        node = node->get_parent();
    }
}
/*---------------------------------------------------------------------------*/
std::shared_ptr<IAction> MCTSAgent::get_random_action(std::shared_ptr<IState> AState) const {
    /** selects an action among the untried ones */
    auto actions = AState->get_actions();
    //randomly pick an action
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, actions.size()-1);
    return actions[distrib(gen)];
}