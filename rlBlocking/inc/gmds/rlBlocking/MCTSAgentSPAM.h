/*---------------------------------------------------------------------------*/
#ifndef MATCHING_MCTSAGENT_H
#define MATCHING_MCTSAGENT_H
/*---------------------------------------------------------------------------*/
#include <memory>
#include <string>
/*---------------------------------------------------------------------------*/
#include "gmds/rlBlocking/IStateSPAM.h"
#include "gmds/rlBlocking/IActionSPAM.h"
#include "gmds/rlBlocking/IRewardFunctionSPAM.h"
#include "gmds/rlBlocking/MCTSTreeSPAM.h"
/*---------------------------------------------------------------------------*/
class MCTSAgent {
 public:

	/**@brief Enumerate type for handling output files in debug mode.
	 */
	enum DEBUG_OUTPUT_MODE{
		OUT_NO,         // no output
		OUT_ITERATION,  // output at some iterations
		OUT_END_ONLY    // output at the end only
	};

	/**@brief Constructor
     *
     * @param ARewardFunction reward function to evaluate a state
     * @param AMaxIterations  max number of iterations
     * @param AMaxSeconds     max number seconds for the process
     * @param AMaxSimulationDepth max depth during the simulation stage
	 */
	MCTSAgent(const IRewardFunction* ARewardFunction,
	          const int AMaxIterations=100000,
	          const int AMaxSeconds=600,
	          const int AMaxSimulationDepth=100);
	/**@brief default destructor
	 */
	virtual ~MCTSAgent();

	/**@brief Activates the export of output files during the process
     *
     * @param AFileNamePrefix file name prefix to output files
     * @param AOutputMode output mode
     * @param AFrequency  output frequency only used for the OUT_ITERATION mode. A frequency of one indicates
     *                    that it export the tree for every mcts cycle (selection, expansion, simulation, backpropagation)
	 */
	void activate_debug_mode(const std::string& AFileNamePrefix= "mcts",
	                         const DEBUG_OUTPUT_MODE AOutputMode = OUT_END_ONLY,
	                         const int AFrequency= 1);

	/**@brief desactivate the debug output mode*/
	void desactivate_debug_output();

	/**@brief Launch the algorithm
     * @param ARootState the state we start from to build the MCSTree
	 */
	void run(std::shared_ptr<IState> ARootState);

	/**@brief Compute the best solution by traversing the tree from its root to a leaf. We consider
     * here the best solution as being the most visited node at each level.
     *
     * @return the "best visited" solution
	 */
	std::shared_ptr<IState> get_best_solution_visited();

	/**@brief Compute the best solution by traversing the tree from its root to a leaf. We consider
     * here the best solution as being the node with the best uct score at each level.
     *
     * @return the "best uct" solution
	 */
	std::shared_ptr<IState> get_best_solution_uct();

	/**@brief Compute the best solution by traversing the tree from its root to a leaf. We consider
     * here the best solution as being the best winrate node at each level.
     *
     * @return the "best winrate" solution
	 */
	std::shared_ptr<IState> get_best_solution_winrate();

	/**@brief Compute the best node. We consider
     * here the best node as being the best uct score.
     *
     * @return the "best node" uct
	 */
	MCTSTree* get_best_node_uct();

	/**@brief provides the number of iterations done by the algorithm
	 */
	int get_nb_iterations() const {return m_nb_iterations;}
	/**@brief provides the number of seconds used to run the algorithm
	 */
	double get_nb_seconds() const {return m_nb_seconds;}

 private:
	/**@brief Among all the possible action generated from @p AState, pick one randomly
    * @param[in] AState stage we generate an action from
    * @return the generated action
	 */
	std::shared_ptr<IAction> get_random_action(std::shared_ptr<IState> AState) const;
 private:

	/**@brief Selection induces a decision policy, known as the tree policy, to navigate
     * through the existing decision tree, attempting to strike a balance between
     * exploration of unknown decision paths and exploitation of known, promising decision paths.
     * @param[in] ANode the node we start the selection from
     * @return the selected node
	 */
	MCTSTree* select(MCTSTree* ANode);
	/**@brief Simulation stage. This process returns an expected reward that is computed using a random
     * strategy to compute next stages.
     *
     * @param[in] ANode the node we simulate from
     * @return the reward obtained by the simulation process
	 */
	double simulate(MCTSTree* ANode);
	/**@brief Creates a new child for @p ANode by applying an untried actions. If there is no action to perform
     * it returns the input node @p ANode
     *
     * @param[in] ANode the node to expand from
     * @return the new child node
	 */
	MCTSTree* expand(MCTSTree* ANode);
	/**@brief Back propagate the reward obtained at node @p ANode to all the nodes met in his parenthood
     *
     * @param[in] ANode the node we start from
     * @param[in] AReward the reward to give to all parents
	 */
	void back_propagate(MCTSTree* ANode, double AReward);

	void export_tree();
 private:
	/** the tree we build during the run() process */
	MCTSTree* m_tree;
	/** the reward function we use to evaluate each state */
	const IRewardFunction* m_reward_function;
	/** the max number of MCTS iterations that are performed*/
	const int m_max_iterations;
	/** the max number of seconds we accept to spend in the process*/
	const int m_max_seconds;
	/** the max size of a simulation */
	int m_simulation_depth;
	/** the number of iterations that are really done*/
	int m_nb_iterations;
	/** the number of seconds really spent in the process*/
	double m_nb_seconds;
	/** debug output activate (true) or not (false) */
	bool m_debug_activate;
	/** debug output prefix file name */
	std::string m_debug_file_prefix;
	/** debug mode*/
	DEBUG_OUTPUT_MODE m_debug_mode;
	/** frequency in tuned mode*/
	int m_debug_frequency;

};





/*---------------------------------------------------------------------------*/
#endif //MATCHING_MCTSAGENT_H
/*---------------------------------------------------------------------------*/
