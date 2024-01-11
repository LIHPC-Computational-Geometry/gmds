/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSTree.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
/*-------------------------------- MCTS NODE ---------------------------------*/
/*----------------------------------------------------------------------------*/
MCTSNode::MCTSNode(gmds::MCTSNode *AParent, const gmds::MCTSMove *AMove, MCTSState *AState)
	:parent(AParent), move(AMove),state(AState), score(0.0), number_of_simulations(0), size(0)
{
	children = new std::vector<MCTSNode *>();
	untried_actions = state->actions_to_try();
	terminal = state->is_terminal();



}
/*----------------------------------------------------------------------------*/
MCTSNode::~MCTSNode() {
	delete state;
	delete move;
	for (auto *child : *children) {
		delete child;
	}
	delete children;
	while (!untried_actions->empty()) {
		delete untried_actions->front();    // if a move is here then it is not a part of a child node and needs to be deleted here
		untried_actions->pop();
	}
	delete untried_actions;
}

/*----------------------------------------------------------------------------*/
void MCTSNode::expand() {
	if (is_terminal()) {              // can legitimately happen in end-game situations
		rollout();                    // keep rolling out, eventually causing UCT to pick another node to expand due to exploration
		return;
	} else if (is_fully_expanded()) {
		std::cerr << "Warning: Cannot expanded this node any more!" << std::endl;
		return;
	}
	// get next untried action
	MCTSMove *next_move = untried_actions->front();     // get value
	untried_actions->pop();                              // remove it
	MCTSState *next_state = state->next_state(next_move);

	if(state->get_quality() == next_state->get_quality()){
		const unsigned int nb_same_quality = this->nb_same_quality + 1;
	}
	else{
		const unsigned int nb_same_quality = 0;
	}
	// build a new MCTS node from it
	MCTSNode *new_node = new MCTSNode(this,next_move,next_state);
	// rollout, updating its stats
	new_node->rollout();
	// add new node to tree
	children->push_back(new_node);
}
/*----------------------------------------------------------------------------*/
const MCTSState *MCTSNode::get_current_state() const
{
	return state;
}
/*----------------------------------------------------------------------------*/
bool
MCTSNode::is_terminal() const
{
	return terminal;
}

/*----------------------------------------------------------------------------*/
const MCTSMove *MCTSNode::get_move() const {
	return move;
}
/*----------------------------------------------------------------------------*/
bool
MCTSNode::is_fully_expanded() const
{
	return is_terminal() || untried_actions->empty();
}
/*----------------------------------------------------------------------------*/
unsigned int MCTSNode::get_size() const {
	return size;
}
/*----------------------------------------------------------------------------*/
MCTSNode *MCTSNode::select_best_child(double c) const {
	/** selects best child based on the winrate of whose turn it is to play */
	if (children->empty()) return NULL;
	else if (children->size() == 1) return children->at(0);
	else {
		double uct, max = -1;
		MCTSNode *argmax = NULL;
		for (auto *child : *children) {
			double winrate = child->score / ((double) child->number_of_simulations);
			if (c > 0) {
				uct = winrate +
				      c * sqrt(log((double) this->number_of_simulations) / ((double) child->number_of_simulations));
			} else {
				uct = winrate;
			}
			if (uct > max) {
				max = uct;
				argmax = child;
			}
		}
		return argmax;
	}
}
/*----------------------------------------------------------------------------*/
void
MCTSNode::rollout()
{
	double w = state->state_rollout();
	backpropagate(w, 1);
}
/*----------------------------------------------------------------------------*/
void MCTSNode::backpropagate(double w, int n) {
	score += w;
	number_of_simulations += n;
	if (parent != NULL) {
		parent->size++;
		parent->backpropagate(w, n);
	}
}

/*----------------------------------------------------------------------------*/
MCTSNode *MCTSNode::advance_tree(const MCTSMove *m) {
	//TODO
	// Find child with this m and delete all others
	MCTSNode *next = NULL;
	for (auto *child: *children) {
		if (*(child->move) == *(m)) {
			next = child;
		} else {
			delete child;
		}
	}
	// remove children from queue so that they won't be re-deleted by the destructor when this node dies (!)
	this->children->clear();
	// if not found then we have to create a new node
	if (next == NULL) {
		// Note: UCT may lead to not fully explored tree even for short-term children due to terminal nodes being chosen
		std::cout << "INFO: Didn't find child node. Had to start over." << std::endl;
		MCTSState *next_state = state->next_state(m);
		next = new MCTSNode(this, m,next_state);
	} else {
		next->parent = NULL;     // make parent NULL
		                         // IMPORTANT: m and next->move can be the same here if we pass the move from select_best_child()
		                         // (which is what we will typically be doing). If not then it's the caller's responsibility to delete m (!)
	}
	// return the next root
	return next;
}
/*----------------------------------------------------------------------------*/
void MCTSNode::print_stats() const {
#define TOPK 10
	if (number_of_simulations == 0) {
		std::cout << "Tree not expanded yet" << std::endl;
		return;
	}
	std::cout << "___ INFO _______________________" << std::endl
	     << "Tree size: " << size << std::endl
	     << "Number of simulations: " << number_of_simulations << std::endl
	     << "Branching factor at root: " << children->size() << std::endl;
	// print TOPK of them along with their winrates
//	std::cout << "Best moves:" << std::endl;
//	for (int i = 0 ; i < children->size() && i < TOPK ; i++) {
//		std::cout << "  " << i + 1 << ". " << children->at(i)->move->sprint() << "  -->  "
//		     << std::setprecision(4) << 100.0 * children->at(i)->calculate_winrate(state->player1_turn()) << "%" << endl;
//	}
	std::cout << "________________________________" << std::endl;
}


/*----------------------------------------------------------------------------*/
double MCTSNode::calculate_winrate() const {
	return score / number_of_simulations;
}

/*----------------------------------------------------------------------------*/
/*--------------------------------  MCTS TREE --------------------------------*/
/*----------------------------------------------------------------------------*/
MCTSTree::MCTSTree(MCTSState* starting_state)
{
	assert(starting_state != NULL);
	root = new MCTSNode(NULL, NULL, starting_state);
}
/*----------------------------------------------------------------------------*/
MCTSTree::~MCTSTree()
{
	delete root;
}
/*----------------------------------------------------------------------------*/
MCTSNode *MCTSTree::select(double c) {
	MCTSNode *node = root;
	while (!node->is_terminal()) {
		if (!node->is_fully_expanded()) {
			return node;
		} else {
			node = node->select_best_child(c);
		}
	}
	return node;
}
/*----------------------------------------------------------------------------*/
void MCTSTree::grow_tree(int max_iter, double max_time_in_seconds) {
	MCTSNode *node;
	double dt;
#ifdef DEBUG
	std::cout << "Growing tree..." << std::endl;
#endif
	time_t start_t, now_t;
	time(&start_t);
	for (int i = 0 ; i < max_iter ; i++){
		// select node to expand according to tree policy
		node = select();
		// expand it (this will perform a rollout and backpropagate the results)
		node->expand();
		// check if we need to stop
		time(&now_t);
		dt = difftime(now_t, start_t);
		if (dt > max_time_in_seconds) {
#ifdef DEBUG
			std::cout << "Early stopping: Made " << (i + 1) << " iterations in " << dt << " seconds." << std::endl;
#endif
			break;
		}
	}
#ifdef DEBUG
	time(&now_t);
	dt = difftime(now_t, start_t);
	cout << "Finished in " << dt << " seconds." << endl;
#endif
}
/*----------------------------------------------------------------------------*/
MCTSNode *MCTSTree::select_best_child() {
	return root->select_best_child(0.0);
}
/*----------------------------------------------------------------------------*/
void MCTSTree::advance_tree(const MCTSMove *move) {
	MCTSNode *old_root = root;
	root = root->advance_tree(move);
	delete old_root;       // this won't delete the new root since we have emptied old_root's children
}

/*----------------------------------------------------------------------------*/
unsigned int MCTSTree::get_size() const {
	return root->get_size();
}
/*----------------------------------------------------------------------------*/
const MCTSState *MCTSTree::get_current_state() const
{
	return root->get_current_state();
}
/*----------------------------------------------------------------------------*/
void MCTSTree::print_stats() const
{
	root->print_stats();
}
/*----------------------------------------------------------------------------*/

