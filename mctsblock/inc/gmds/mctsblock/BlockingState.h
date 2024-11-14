/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_BLOCKING_STATE_H
#define GMDS_MCTS_BLOCKING_STATE_H
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include "mcts/IState.h"
#include <LIB_GMDS_MCTSBLOCK_export.h>
#include <gmds/mctsblock/Blocking.h>
#include <iostream>
#include <queue>
#include <memory>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace mctsblock {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_MCTSBLOCK_API BlockingState : public IState
{
 public:


	static  double weight_nodes;
	static  double weight_edges;
	static  double weight_faces;

	BlockingState(std::shared_ptr<Blocking> AB, int ADepth = 0, std::deque<double> APrevScore = std::deque<double>());
	BlockingState(const BlockingState &AState);

	std::vector<std::shared_ptr<IAction>> get_actions() const override;

	bool is_terminal() const override;

	bool win() const override;

		bool lost() const override;
	   bool  draw() const override;
	std::string write(const std::string &AFileName, const int AStageIndex, const int ANodeId, const int ADepth) const override;


	/**@brief computes the score of the current blocking and maintain
	 * the memory stack of scores
	 *
	 * @return
	 */
	double computeScore();
	void updateMemory(const double AScore);

	/**@brief computes the score of the classification
	 *
	 * @return
	 */
	double computeScoreClassification();
	/**@brief computes the min edge length
	 *
	 * @return
	 */
	double computeMinEdgeLenght() const;


	std::shared_ptr<Blocking> get_blocking() const {return m_blocking;}
	int get_depth() const {return m_depth;}
	/**
	 * Stack of memorized score. The last score in the stack (back) is the score of the current state
	 * @return
	 */
	std::deque<double>  get_memory() const {return m_memory_scores;}

	double get_expected_optimal_score() const {return m_expected_optimal_score;}
 private:

	/**@brief This method return all the possible cut
	 * @return a vector with only pair in, the first (pair.first) is the edge, and the second (pair.second) is the param to cut
	 */
	std::vector<std::shared_ptr<IAction>> get_possible_cuts() const;

	/**@brief This method returns all the possible block erasing. A block can be erased if
	 *  - it doesn't not belong a corner that is the only one to capture a point.
	 *  - its centroid is not inside the volume
	 * @return a vector of block ids. */
	std::vector<std::shared_ptr<IAction>> get_possible_block_removals() const;

 private:
	/** the memory depth we store*/
	static const int m_memory_depth;

	double m_expected_optimal_score;
	/** the blocking we act on. Note that the geometric model is known by the blocking */
	std::shared_ptr<Blocking> m_blocking;

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
