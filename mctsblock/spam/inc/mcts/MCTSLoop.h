/*---------------------------------------------------------------------------*/
#ifndef MATCHING_MCTSLOOP_H
#define MATCHING_MCTSLOOP_H
/*---------------------------------------------------------------------------*/
#include <memory>
#include <iostream>
/*---------------------------------------------------------------------------*/
#include "mcts/IState.h"
#include "mcts/IAction.h"
#include "mcts/IRewardFunction.h"
#include "mcts/MCTSAgent.h"
/*---------------------------------------------------------------------------*/
class MCTSLoop {
public:

    /**@brief Enumerate type for node selection at the end of each loop
     * stage
     */
    enum NODE_SELECTION {
        BEST_SOLUTION, // pick the best child in the full tree
        BEST_CHILD // pick the best direct child of the current node
    };

    /**@brief Constructor
     *
     * @param AAgent the agent we want to train
     * @param ARootState the environment state we work on
     * @param AMaxIterations  max number of iterations
     * @param AMaxSeconds     max number seconds for the whole process
     * @param ADisplayInfo      display info in std::cout stream if true
     */
    MCTSLoop(MCTSAgent& AAgent,
             std::shared_ptr<IState> AState,
             const NODE_SELECTION ASelectionMode = BEST_CHILD,
             const int AMaxIterations=100000,
             const bool ADisplayInfo = false);
    /**@brief default destructor
     */
    virtual ~MCTSLoop();

    /**@brief Launch the loop procedure
     */
    void run();

    /**@brief Gives the number of iterations used during the last run
     */
    int get_nb_iterations() const {return m_nb_iterations;}

    /**@brief Gives access to the used agent */
    MCTSAgent& get_agent() const {return m_agent;}

private:
    /** the agent we use */
    MCTSAgent& m_agent;
    /** Root state we start from */
    std::shared_ptr<IState> m_root_state;
    /** mode selection for picking a node at the end of a loop iteration*/
    NODE_SELECTION m_selection_mode;
    /** the max number of loop iterations that are performed*/
    const int m_max_iterations;
    /** the number of iterations that are really done*/
    int m_nb_iterations;
    /** display info or not on the std::cout stream */
    bool m_display_info;
};

/*---------------------------------------------------------------------------*/
#endif //MATCHING_MCTLOOP_H
/*---------------------------------------------------------------------------*/
