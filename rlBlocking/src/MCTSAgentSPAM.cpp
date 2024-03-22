/*---------------------------------------------------------------------------*/
#include "gmds/rlBlocking/MCTSAgentSPAM.h"
#include <chrono>
#include <iostream>
#include <random>
#include <fstream>
#include <nlohmann/json.hpp>
/*---------------------------------------------------------------------------*/
using json = nlohmann::json;
/*---------------------------------------------------------------------------*/

MCTSAgent::MCTSAgent(const IRewardFunction* ARewardFunction,
                     const int AMaxIterations,
                     const int AMaxSeconds,
                     const int AMaxSimulationDepth)
  : m_tree(nullptr),
  m_reward_function(ARewardFunction),
  m_max_iterations(AMaxIterations),
  m_max_seconds(AMaxSeconds),
  m_simulation_depth(AMaxSimulationDepth),
  m_debug_activate(false),
  m_debug_file_prefix("mcts"),
  m_debug_mode(MCTSAgent::OUT_END_ONLY),
  m_debug_frequency(1)
{}
/*---------------------------------------------------------------------------*/
MCTSAgent::~MCTSAgent() {
	delete m_tree;
}
/*---------------------------------------------------------------------------*/
void MCTSAgent::run(std::shared_ptr<IState> ARootState) {
	//check that an existing tree was not used during a previous iteration
	if(m_tree)
		delete m_tree;
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
		if(m_debug_activate && m_debug_mode==OUT_ITERATION){
			if(i%m_debug_frequency==0)
				export_tree();
		}
	}
	m_nb_iterations=i;
	m_nb_seconds =elapsed.count();
	if (m_debug_activate && (m_debug_mode==OUT_END_ONLY || m_debug_mode==OUT_ITERATION)){
		export_tree();
	}
}
/*---------------------------------------------------------------------------*/
std::shared_ptr<IState> MCTSAgent::get_best_solution_visited() {
	const MCTSTree* node = m_tree;
	while (!node->is_terminal() && node->has_children()){
		node=node->get_most_visited_child();
	}
	return node->get_state();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MCTSTree* MCTSAgent::get_best_node_uct(){
	const MCTSTree* node = m_tree;

	return node->get_best_uct_child();
}
/*---------------------------------------------------------------------------*/
std::shared_ptr<IState> MCTSAgent::get_best_solution_uct(){
	const MCTSTree* node = m_tree;
	while (!node->is_terminal() && node->has_children()){
		node=node->get_best_uct_child();
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

	//if we are here, there is an issue
	throw std::runtime_error("Error: expansion has note been done");
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
/*---------------------------------------------------------------------------*/
void MCTSAgent::activate_debug_mode(const std::string &AFileNamePrefix,
                               const DEBUG_OUTPUT_MODE AOutputMode,
                               const int AFrequency)
{
	m_debug_activate = true;
	m_debug_file_prefix = AFileNamePrefix;
	m_debug_mode = AOutputMode;
	m_debug_frequency=AFrequency;
}
/*---------------------------------------------------------------------------*/
void MCTSAgent::desactivate_debug_output() {
	m_debug_activate=false;
}
/*---------------------------------------------------------------------------*/
void MCTSAgent::export_tree() {
	static int file_index=0;
	json j;
	std::vector<MCTSTree*> to_do;
	to_do.push_back(m_tree);
	while (!to_do.empty()){
		//get the last node
		auto n = to_do.back();
		to_do.pop_back();

		//export the state of n
		n->get_state()->write(m_debug_file_prefix,file_index,
		                      n->get_id(),n->get_depth());

		j["nodes"].push_back(json{{"id",n->get_id()},
		                           {"depth",n->get_depth()},
		                           {"reward",n->cumulative_reward},
		                           {"visits",n->number_visits}});

		if(n->get_parent()!= nullptr){
			//means n is not the root
			j["links"].push_back(json{{"parent",n->get_parent()->get_id()},
			                           {"child",n->get_id()}});

		}
		auto children = n->get_children();
		to_do.insert(to_do.end(),children.begin(),children.end());
	}
	std::ofstream file;
	file.open (m_debug_file_prefix+"_"+std::to_string(file_index)+".json");
	file<<j;
	file.close();
	file_index++;
}