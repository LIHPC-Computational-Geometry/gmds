/*---------------------------------------------------------------------------*/
#include <cmath>
#include <iostream>
#include <cfloat>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <memory>
/*---------------------------------------------------------------------------*/
#include "mcts/MCTSTree.h"
/*---------------------------------------------------------------------------*/
MCTSTree::MCTSTree(std::shared_ptr<IState> AState,
                   std::shared_ptr<IAction> AAction,
                   MCTSTree *AParent)
: m_state(AState), m_action(AAction), m_parent(AParent),
  cumulative_reward(0), number_visits(0), number_win(0),
  number_draw(0), number_lost(0), sq_cumulative_reward(0)
{
    static int global_id=0;
    m_id      = global_id++;
    m_depth   = (m_parent==nullptr)?0:m_parent->m_depth+1;
    m_actions = m_state->get_actions();
}
/*---------------------------------------------------------------------------*/
MCTSTree::~MCTSTree() {
    for(auto c: m_children)
        delete c;
}
/*---------------------------------------------------------------------------*/
int MCTSTree::get_id() const {return m_id;}
/*---------------------------------------------------------------------------*/
int MCTSTree::get_depth() const {return m_depth;}
/*---------------------------------------------------------------------------*/
std::shared_ptr<IState> MCTSTree::get_state() const {
    return m_state;
}
/*---------------------------------------------------------------------------*/
std::shared_ptr<IAction> MCTSTree::get_action() const {
    return m_action;
}
/*---------------------------------------------------------------------------*/
MCTSTree* MCTSTree::get_parent() const {
    return  m_parent;
}
/*---------------------------------------------------------------------------*/
bool MCTSTree::is_terminal() const {
    return m_state->is_terminal();
}
/*---------------------------------------------------------------------------*/
bool MCTSTree::is_fully_expanded() const {
    return m_children.empty() == false &&
           m_children.size() == m_actions.size();
}

/*---------------------------------------------------------------------------*/
MCTSTree* MCTSTree::get_most_visited_child() const {
    int most_visits = -1;
    MCTSTree* best_node = nullptr;

    // iterate all  children and find most visited
    for(auto c: m_children) {
		  std::cout<<"Child "<<c<<" from action "<<c->get_action()->get_description()<<" (visits, win, lost, draw, cumul. reward): "<<c->number_visits<<", "<<c->number_win<<", "<<c->number_lost<<", "<<c->number_draw<<", "<<c->cumulative_reward<<std::endl;

        if(c->number_visits > most_visits) {
            most_visits = c->number_visits;
            best_node = c;
        }
    }
    if(best_node== nullptr)
        throw std::runtime_error("Error when visiting children");
	 std::cout<<"--> Best "<<best_node<<" (visits, win, lost, draw, cumul. reward): "
	           <<best_node->number_visits<<", "
	           <<best_node->number_win<<", "
	           <<best_node->number_lost<<", "
	           <<best_node->number_draw<<", "
	           <<best_node->cumulative_reward<<std::endl;

    return best_node;
}
/*---------------------------------------------------------------------------*/
MCTSTree* MCTSTree::get_child(const int AI) const {
    return m_children[AI];
}

std::vector<MCTSTree *> MCTSTree::get_children() const {
    return m_children;
}
/*---------------------------------------------------------------------------*/
bool MCTSTree::has_children() const {
    return !m_children.empty();
}
/*---------------------------------------------------------------------------*/
MCTSTree*  MCTSTree::expand()  {
    // sanity check that we're not already fully expanded
    if(is_fully_expanded())
        throw std::runtime_error("Error: cannot expand a fully expanded node");

    if (is_terminal())
        return this;

    // add the next action in queue as a child
    return add_child_with_action(m_actions[m_children.size()]);

}

//--------------------------------------------------------------
// create a clone of the current state, apply action, and add as child
MCTSTree*  MCTSTree::add_child_with_action(std::shared_ptr<IAction> AAction) {
    // create a new TreeNode with the same state (will get cloned) as this TreeNode
    auto next_state = AAction->apply_on(m_state);

    MCTSTree* child_node = new MCTSTree(next_state, AAction, this);

    m_children.push_back(child_node);

    return child_node;

}