/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSTREE_H
#define GMDS_MCTSTREE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSState.h>
#include <gmds/rlBlocking/MCTSMove.h>
#include <iostream>
#include <cmath>
#include <cassert>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSNode
 *  @brief  Class that provides ....
 */
class LIB_GMDS_RLBLOCKING_API MCTSNode {
	/** @brief Yes if the node have no childs */
	bool terminal;
	/** @brief Number of nodes in the tree from the node.  */
	unsigned int size;
	/** @brief number of parent nodes having the same quality*/
	unsigned int nb_same_quality;
	/** @brief Number of visits*/
	unsigned int number_of_simulations;
	/** @brief  e.g. number of wins (could be int but double is more general if we use evaluation functions)*/
	double score;
	/** @brief state of the current node */
	MCTSState *state;
	/** @brief  move to get here from parent node's state*/
	const MCTSMove *move;
	/** @brief the chilren for the current node */
	std::vector<MCTSNode *> *children;
	/** @brief  the parent for the current node*/
	MCTSNode *parent;
	/** @brief  queue of untried actions*/
	std::deque<MCTSMove *> *untried_actions;
	/** @brief  update the nb simulations and the score after a rollout*/
	void backpropagate(double w, int n);
 public:

	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
         * @param AParent the parent the node
         * @param	AMove the action to access at this node
	 */
	MCTSNode(MCTSNode *AParent, const MCTSMove *AMove, MCTSState *AState);

	/*------------------------------------------------------------------------*/
	/** @brief  Destructor.	*/
	virtual ~MCTSNode();

	/** @brief  Check if the node is fully expanded	*/
	bool is_fully_expanded() const;
	/** @brief  Check if the node is terminal.
	 * @param	Number max of parents nodes with the same quality
	 * */
	bool is_terminal() const;
	/** @brief  Return the different moves/actions possible for a node	*/
	const MCTSMove *get_move() const;
	/** @brief  Return the size.	*/
	unsigned int get_size() const;
	/** @brief  Expand the node.	*/
	void expand();
	/** @brief  Do a rollout.	*/
	void rollout();
	/** @brief  Select the most promising child of the root node	*/
	MCTSNode *select_best_child(double c) const;
	/** @brief  Find child with this m and delete all others.
	 * @param m the selected move
	 * @return the next root*/
	MCTSNode *advance_tree(const MCTSMove *m);
	/** @brief  Return the state of the node.	*/
	const MCTSState *get_current_state() const;
	/** @brief  Return the children of the node.	*/
	std::vector<MCTSNode *> *get_children();
	/** @brief  Print the tree and the stats.	*/
	void print_stats() const;
	/** @brief  Calculate the q rate of a node. It's: wins-looses	*/
	double q_rate() const;
	/** @brief  Calculate UCT. 	*/
	double calculate_UCT() const;
	/** @brief  Calculate winrate. 	*/
	double calculate_winrate() const;



 private:
	/** a mesh */
	//Mesh* m_mesh;
};

/*----------------------------------------------------------------------------------------*/
/** @class  MCTSTree
 *  @brief  Class that provides ....
 */
class LIB_GMDS_RLBLOCKING_API MCTSTree
{
	MCTSNode *root;
 public:

	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
         * @param
	 */
	MCTSTree(MCTSState *starting_state);

	/*------------------------------------------------------------------------*/
	/** @brief  Destructor.	*/
	virtual ~MCTSTree();

	/** @brief select child node to expand according to tree policy (UCT).
         * @param c  exploration parameter, theoretically equal to √2
         * @return a node
	 */
	MCTSNode *select(double c=1.41);
	MCTSNode *select_best_child();
	void grow_tree(int max_iter, double max_time_in_seconds);
	/** @brief if the move is applicable advance the tree, else start over
         * @param move the move to do
         * */
	void advance_tree(const MCTSMove *move);
	/** @brief get the size of the tree. */
	unsigned int get_size() const;
	/** @brief get the size of the tree. */
	const MCTSState *get_current_state() const;
	/** @brief Print stats. */
	void print_stats() const;


 private:
	/** a mesh */
	//Mesh* m_mesh;
};
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSTREE_H
/*----------------------------------------------------------------------------------------*/
