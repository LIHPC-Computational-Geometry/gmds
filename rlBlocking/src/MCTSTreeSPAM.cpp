/*---------------------------------------------------------------------------*/
#include <cmath>
#include <iostream>
#include <cfloat>
#include <limits>
#include <random>
#include <algorithm>
#include <memory>
/*---------------------------------------------------------------------------*/
#include "gmds/rlBlocking/MCTSTreeSPAM.h"
/*---------------------------------------------------------------------------*/
MCTSTree::MCTSTree(std::shared_ptr<IState> AState, std::shared_ptr<IAction> AAction, MCTSTree *AParent)
: m_state(AState), m_action(AAction), m_parent(AParent),
  cumulative_reward(0), number_visits(0)
{
    m_actions =m_state->get_actions();
}
/*---------------------------------------------------------------------------*/
MCTSTree::~MCTSTree() {
    for(auto c: m_children)
        delete c;
}
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
        if(c->number_visits > most_visits) {
            most_visits = c->number_visits;
            best_node = c;
        }
    }
    return best_node;



}
/*---------------------------------------------------------------------------*/
MCTSTree* MCTSTree::get_child(const int AI) const {
    return m_children[AI];
}
/*---------------------------------------------------------------------------*/
bool MCTSTree::has_children() const {
    return !m_children.empty();
}
/*---------------------------------------------------------------------------*/
MCTSTree*  MCTSTree::expand()  {
    // sanity check that we're not already fully expanded
    if(is_fully_expanded())
        return nullptr;

    if (is_terminal())
        return this;

    // if this is the first expansion and we haven't yet got all of the possible actions
    if(m_actions.empty()) {
        // retrieve list of actions from the state
        m_actions = m_state->get_actions();

        // randomize the order
        std::shuffle(m_actions.begin(),
                     m_actions.end(),
                     std::mt19937(std::random_device()()));
    }

    // add the next action in queue as a child
    return add_child_with_action(m_actions[m_children.size()]);

}

//--------------------------------------------------------------
// create a clone of the current state, apply action, and add as child
MCTSTree*  MCTSTree::add_child_with_action(std::shared_ptr<IAction> AAction) {
    // create a new TreeNode with the same state (will get cloned) as this TreeNode
    auto next_state = m_state->apply(AAction);

    MCTSTree* child_node = new MCTSTree(next_state, AAction, this);

    m_children.push_back(child_node);

    return child_node;

}
/*---------------------------------------------------------------------------*/
MCTSTree* MCTSTree::get_best_uct_child(double AC) const {
    // sanity check
    if(!is_fully_expanded())
        return nullptr;

    float best_utc_score = -std::numeric_limits<float>::max();
    MCTSTree* best_node = nullptr;

    // iterate all immediate children and find best UTC score
    int num_children = m_children.size();
    for(auto i = 0; i < num_children; i++) {
        auto child = m_children[i];
        float uct_exploitation = (float)child->cumulative_reward / (child->number_visits + FLT_EPSILON);
        float uct_exploration = sqrt(log((float)this->number_visits + 1) / (child->number_visits + FLT_EPSILON) );
        float uct_score = uct_exploitation + AC * uct_exploration;

        if(uct_score > best_utc_score) {
            best_utc_score = uct_score;
            best_node = child;
        }
    }

    return best_node;

}