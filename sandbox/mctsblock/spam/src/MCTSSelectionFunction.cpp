/*---------------------------------------------------------------------------*/
#include <cmath>
#include <limits>
#include <stdexcept>
/*---------------------------------------------------------------------------*/
#include "mcts/MCTSTree.h"
#include "mcts/MCTSSelectionFunction.h"
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
		      auto np = (double)ANode->number_visits + 1;
		      auto ni = (double)child->number_visits + FLT_EPSILON;
		      auto vi =  (double)child->cumulative_reward;
		      auto utc_score= vi/ni + m_c * sqrt(2*log(np) / ni);

		      if(utc_score > best_utc_score) {
			      best_utc_score = utc_score;
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
		      auto np = (double)ANode->number_visits + 1;
		      auto ni = (double)child->number_visits + FLT_EPSILON;
		      auto vi =  (double)child->cumulative_reward;
		      auto utc= vi/ni + m_c * sqrt(2*log(np) / ni);

		      auto sum_x2 = child->sq_cumulative_reward;
		      auto utc_single = sqrt((sum_x2 - ni* pow(vi/ni,2) +m_d) / ni);
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
