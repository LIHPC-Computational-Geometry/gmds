/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSAgent.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
MCTSAgent::MCTSAgent(gmds::MCTSState *starting_state, int max_iter, int max_seconds, int max_same_quality)
	:max_iter(max_iter),max_seconds(max_seconds),max_same_quality(max_same_quality)
{
	tree = new MCTSTree(starting_state);
}

/*----------------------------------------------------------------------------*/
MCTSAgent::~MCTSAgent(){
	delete tree;
}
/*----------------------------------------------------------------------------*/
const MCTSMove *MCTSAgent::genmove()
{
	// If game ended from opponent move, we can't do anything
	if (tree->get_current_state()->is_terminal()) {
		return NULL;
	}
#ifdef DEBUG
	std::cout << "___ DEBUG ______________________" << endl
	     << "Growing tree..." << std::endl;
#endif
	tree->grow_tree(max_iter, max_seconds);
#ifdef DEBUG
	cout << "Tree size: " << tree->get_size() << endl
	     << "________________________________" << endl;
#endif
	MCTSNode *best_child = tree->select_best_child();
	if (best_child == NULL) {
		std::cerr << "Warning: Tree root has no children! Possibly terminal node!" << std::endl;
		return NULL;
	}
	const MCTSMove *best_move = best_child->get_move();
	tree->advance_tree(best_move);
	return best_move;
}
/*----------------------------------------------------------------------------*/
const
   MCTSState *MCTSAgent::get_current_state() const
{
	return tree->get_current_state();
}
/*----------------------------------------------------------------------------*/
