#ifndef GMDS_MCTSAGENT_H
#define GMDS_MCTSAGENT_H

#include <gmds/rlBlocking/MCTSTree.h>
#include <iostream>
#include <random>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
class LIB_GMDS_RLBLOCKING_API MCTSAgent
{
	// example of an agent based on the MCTS_tree. One can also use the tree directly.
	MCTSTree *tree;
	int max_iter, max_seconds, max_same_quality;
	std::string m_NameGeo;

 public:
	MCTSAgent(MCTSState *starting_state,std::string ANameGeo, int max_iter = 100000, int max_seconds = 30, int max_same_quality=3);
	~MCTSAgent();
	const MCTSMove *genmove();
	const MCTSState *get_current_state() const;
	void feedback() const {tree->print_stats();}
};
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSAGENT_H
/*----------------------------------------------------------------------------------------*/
