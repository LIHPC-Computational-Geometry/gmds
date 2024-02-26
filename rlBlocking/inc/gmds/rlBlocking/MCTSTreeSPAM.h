/*---------------------------------------------------------------------------*/
#ifndef MATCHING_MCTSTREE_H
#define MATCHING_MCTSTREE_H
/*---------------------------------------------------------------------------*/
#include <memory>
#include "gmds/rlBlocking/IStateSPAM.h"
#include "gmds/rlBlocking/IActionSPAM.h"
#include "gmds/rlBlocking/IRewardFunctionSPAM.h"
/*---------------------------------------------------------------------------*/
class MCTSTree{
public:
    MCTSTree(std::shared_ptr<IState> AState, std::shared_ptr<IAction> AAction = nullptr,
             MCTSTree* AParent = nullptr);

    /**@brief Default destructor
     */
    ~MCTSTree();

    /**@brief Provides the state associated to the root node
     * @return the current root state
     */
    std::shared_ptr<IState> get_state() const;

    /**@brief Give access to the parent node
     * @return Returns the parent node. Returns nullptr for the root node
     */
   MCTSTree* get_parent() const;

    /**@brief Provides the action that allows to create this node from its
     *        parent
     * @return an action, which can be null for the root
     */
    std::shared_ptr<IAction> get_action() const;

    /** Return the child that can be reached applying the action @p AAction.
     * It returns the node if it exists, but do not create it otherwise.
     * @param[in] AAction
     * @return the expected node, or nullptr if it doesn't exist
     */
    MCTSTree* get_child(const int AI) const;

    MCTSTree* get_most_visited_child() const;

    MCTSTree*  add_child_with_action(std::shared_ptr<IAction> new_action);

    MCTSTree*  expand();
    /** @return true if it doesn't remain actions to perform for generating child,
     * false otherwise
     * All children have been expanded and simulated
     */
    bool is_fully_expanded() const;
    /**@brief a node is terminal when the associated state is said "win" or "lost".
     * In other word, the game is over
     * @return true if the game is over, false otherwise.
     */
    bool is_terminal() const;

    bool has_children() const;
    /**@brief Among all the children, provide the action that reach it
     *        based on the exploration/exploitation principle (UCT criterion here)
       * @return the selected action
       */
    MCTSTree*  get_best_uct_child(double AC=1.41) const;


    /** total reward of the current node */
    double cumulative_reward;
    /** number of times the node has been explored */
    int number_visits;

private:
    /** current state **/
    std::shared_ptr<IState> m_state;
    /** action that generates this node from the parent **/
    std::shared_ptr<IAction> m_action;
    /** parent node that is shared with all the siblings**/
    MCTSTree*  m_parent;
    /** children nodes, that are potentially empty. When a children node
     * is generated, it means that an action is remove from the m_untried_actions
     * collection **/
    std::vector<MCTSTree*> m_children;

    std::vector<std::shared_ptr<IAction>> m_actions;
};
/*---------------------------------------------------------------------------*/
#endif //MATCHING_MCTSTREE_H
/*---------------------------------------------------------------------------*/
