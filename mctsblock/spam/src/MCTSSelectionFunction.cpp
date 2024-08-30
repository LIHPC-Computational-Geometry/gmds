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
            float uct_exploitation = (float)child->number_win / (child->number_visits + FLT_EPSILON);
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
SPUCTSelectionFunction::SPUCTSelectionFunction(const double AC, const double AD)
  :m_c(AC), m_d(AD){;}
/*---------------------------------------------------------------------------*/
MCTSTree* SPUCTSelectionFunction::select(MCTSTree *ANode)  const  {
	     // sanity check
	     if(!ANode->is_fully_expanded())
		      throw std::runtime_error("Error: cannot compute the best child for a partially expanded node");

	     float best_score = -std::numeric_limits<float>::max();
	     MCTSTree* best_node = nullptr;

	     // iterate all immediate children and find best UTC score
	     auto children = ANode->get_children();
	     auto num_children = children.size();

	     for(auto i = 0; i < num_children; i++) {
		      auto child = children[i];
		      auto tN = (double)ANode->number_visits + 1;
		      auto tNi = (double)child->number_visits + FLT_EPSILON;
		      auto w =  (double)child->number_win;
		      auto utc= w/tNi + m_c * sqrt(log(tN) / tNi);

		      auto sum_x2 = child->sq_cumulative_reward;
		      auto utc_single = sqrt((sum_x2 - tNi* pow(w/tNi,2) +m_d) / tNi);
		      auto utc_score = utc + utc_single;

		      if(utc_score > best_score) {
			       best_score = utc_score;
			       best_node = child;
		      }
	     }
	     if(best_node== nullptr)
		      throw std::runtime_error("Error when getting the best child of a node");

	     return best_node;
};
/*---------------------------------------------------------------------------*/
