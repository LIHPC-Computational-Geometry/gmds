/*---------------------------------------------------------------------------*/
#ifndef MATCHING_MCTSAGENT_H
#define MATCHING_MCTSAGENT_H
/*---------------------------------------------------------------------------*/
#include <memory>
/*---------------------------------------------------------------------------*/
#include "gmds/rlBlocking/IStateSPAM.h"
#include "gmds/rlBlocking/IActionSPAM.h"
#include "gmds/rlBlocking/IRewardFunctionSPAM.h"
#include "gmds/rlBlocking/MCTSTreeSPAM.h"
/*---------------------------------------------------------------------------*/
class MCTSAgent {
public:

    MCTSAgent(const IRewardFunction* ARewardFunction,
              const int AMaxIterations=100000, const int AMaxSeconds=600,
              const int AMaxSimulationDepth=100);
    virtual ~MCTSAgent();
    void run(std::shared_ptr<IState> ARootState);
    std::shared_ptr<IState> get_best_solution();
    int get_nb_iterations() const {return m_nb_iterations;}
    double get_nb_seconds() const {return m_nb_seconds;}

    /**@brief Among all the possible action generated from
    * @p AState, pick one randomly
    * @param[in] AState stage we generate an action from
    * @return the generated action
    */
    std::shared_ptr<IAction> get_random_action(std::shared_ptr<IState> AState) const;

private:

    /** Selection induces a decision policy, known as the tree policy, to navigate
     * through the existing decision tree, attempting to strike a balance between
     * exploration of unknown decision paths and exploitation of known, promising decision paths.
     * @param[in] ANode the node we start the selectio from
     * @return the selected node and the reward associated to this node
     */
    MCTSTree* select(MCTSTree* ANode);
    double simulate(MCTSTree* ANode);
    MCTSTree* expand(MCTSTree* ANode);
    void back_propagate(MCTSTree* ANode, double AReward);
private:
    MCTSTree* m_tree;
    const IRewardFunction* m_reward_function;
    const int m_max_iterations;
    const int m_max_seconds;
    int m_nb_iterations;
    double m_nb_seconds;
    int m_simulation_depth;
};




/*---------------------------------------------------------------------------*/
#endif //MATCHING_MCTSAGENT_H
/*---------------------------------------------------------------------------*/
