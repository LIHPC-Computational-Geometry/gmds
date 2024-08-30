/*---------------------------------------------------------------------------*/
#ifndef MATCHING_MCTSTREE_H
#define MATCHING_MCTSTREE_H
/*---------------------------------------------------------------------------*/
#include <memory>
#include "mcts/IState.h"
#include "mcts/IAction.h"
#include "mcts/IRewardFunction.h"
/*---------------------------------------------------------------------------*/
class MCTSTree{
public:
    MCTSTree(std::shared_ptr<IState> AState,
             std::shared_ptr<IAction> AAction = nullptr,
             MCTSTree* AParent = nullptr);

    /**@brief Default destructor
     */
    ~MCTSTree();
    /**
     * @return root node id
     */
    int get_id() const;
    /**
     * @return root node depth
     */
    int get_depth() const;
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
    /** Return access to the container of children of (*this). Carefully a
     * handle this access.
     */
    std::vector<MCTSTree*> get_children() const;

    /**@brief Returns the most visited child of the current root node
     * @return a node, or subtree
     */
    MCTSTree* get_most_visited_child() const;

    /**@brief Create a new child that will be obtained by applying an untried action
    * @return a child  node obtained from applying an untried action on the current state
    */
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
    /**@brief Check if the node has children
     * @return true if the node has at least one child, false otherwise
     */
    bool has_children() const;

    /** total reward of the current node */
	 double cumulative_reward;
	 double sq_cumulative_reward;
    /** number of times the node has been explored */
    int number_visits;
	 /** number of times the has has lead to a win */
	 int number_win;
	 /** number of times the has has lead to a draw */
	 int number_draw;
	 /** number of times the has has lead to a lost */
	 int number_lost;

private:
    /**@brief add a child that will be obtained by applying action @p AAction. We do not check if this action
     * has already be performed. The current tree implementation ensures that it is the case.
     *
     * @param[in] AAction the action to perform
     * @return the node obtained from applying @p AAaction on the current state
     */
    MCTSTree*  add_child_with_action(std::shared_ptr<IAction> AAction);

private:
    /** node id that is generated from the parent id*/
    int m_id;
    /** depth in the main tree, stored and note computed. A depth of 0 means tree root*/
    int m_depth;
    /** current state **/
    std::shared_ptr<IState> m_state;
    /** action that generates this node from the parent **/
    std::shared_ptr<IAction> m_action;
    /** parent node that is shared with all the siblings**/
    MCTSTree*  m_parent;
    /** children nodes, that are potentially empty. When a children node
     * is generated, it means that an action from m_actions has been applied*/
    std::vector<MCTSTree*> m_children;

    /** actions to be applied. This list is only generated in the constructor **/
    std::vector<std::shared_ptr<IAction>> m_actions;
};
/*---------------------------------------------------------------------------*/
#endif //MATCHING_MCTSTREE_H
/*---------------------------------------------------------------------------*/
