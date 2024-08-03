/*---------------------------------------------------------------------------*/
#include <cfloat>
#include <mcts/MCTSTree.h>
#include <mcts/MCTSSelectionFunction.h>
/*---------------------------------------------------------------------------*/
UCBSelectionFunction::UCBSelectionFunction(const double AC):m_c(AC){;}
/*---------------------------------------------------------------------------*/
MCTSTree* UCBSelectionFunction::select(MCTSTree *ANode)  const  {
        // sanity check
        if(!ANode->is_fully_expanded())
            throw std::runtime_error("Error: cannot compute the best child for a partially expanded node");

        float best_utc_score = -std::numeric_limits<float>::max();
        MCTSTree* best_node = nullptr;

        // iterate all immediate children and find best UTC score
        auto children = ANode->get_children();
        auto num_children = children.size();
        for(auto i = 0; i < num_children; i++) {
            auto child = children[i];
            float uct_exploitation = (float)child->cumulative_reward / (child->number_visits + FLT_EPSILON);
            float uct_exploration = sqrt(log((float)ANode->number_visits + 1) / (child->number_visits + FLT_EPSILON) );
            float uct_score = uct_exploitation + m_c * uct_exploration;

            if(uct_score > best_utc_score) {
                best_utc_score = uct_score;
                best_node = child;
            }
        }
        if(best_node== nullptr)
            throw std::runtime_error("Error when getting the best child of a node");

        return best_node;
};
/*---------------------------------------------------------------------------*/
