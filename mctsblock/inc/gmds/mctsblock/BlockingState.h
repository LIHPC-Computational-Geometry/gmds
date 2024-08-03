/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_BLOCKING_STATE_H
#define GMDS_MCTS_BLOCKING_STATE_H
/*----------------------------------------------------------------------------*/
#include "mcts/IState.h"
#include <LIB_GMDS_MCTSBLOCK_export.h>
#include <gmds/mctsblock/Blocking.h>
#include <iostream>
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace mctsblock {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_MCTSBLOCK_API BlockingState : public IState
{
 public:
	BlockingState(Blocking *AB, int ADepth = 0, std::deque<double> APrevScore = std::deque<double>());
	BlockingState(const BlockingState &AState);

	std::vector<std::shared_ptr<IAction>> get_actions() const override;

	bool is_terminal() const override;

	bool win() const override;

	std::string write(const std::string &AFileName, const int AStageIndex, const int ANodeId, const int ADepth) const override;

	bool lost() const;

	/**@brief computes the score of the current blocking and maintain
	 * the memory stack of scores
	 *
	 * @return
	 */
	double computeScore();

	void updateMemory(const double AScore);

	Blocking* get_blocking() const {return m_blocking;}
	int get_depth() const {return m_depth;}
	std::deque<double>  get_memory() const {return m_memory_scores;}
 private:

	/**@brief This method return all the possible cut
	 * @return a vector with only pair in, the first (pair.first) is the edge, and the second (pair.second) is the param to cut
	 */
	std::vector<std::shared_ptr<IAction>> get_possible_cuts() const;

	/**@brief This method return all the possible block erasing. A block can be erase if it doesn't not belong a corner
	 * that is the only one to capture a point
	 * @return a vector of block ids. */
	std::vector<std::shared_ptr<IAction>> get_possible_block_removals() const;

 private:
	/** the memory depth we store*/
	static const int m_memory_depth;

	static const double m_weight_nodes;
	static const double m_weight_edges;
	static const double m_weight_faces;

	/** the blocking we act on. Note that the geometric model is known by the blocking */
	Blocking *m_blocking;

	std::set<TCellID> m_boundary_node_ids;
	std::set<TCellID> m_boundary_edge_ids;
	std::set<TCellID> m_boundary_face_ids;
	/** we need to know at what depth we are in the tree */
	int m_depth;
	/** we store the score of the current blocking and the m_memory_depth previsious one.*/
	std::deque<double> m_memory_scores;
};
/*----------------------------------------------------------------------------*/
}     // namespace mctsblock
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MCTS_BLOCKING_STATE_H
