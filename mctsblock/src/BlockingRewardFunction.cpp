/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/BlockingRewardFunction.h>
#include <gmds/mctsblock/BlockingState.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::mctsblock;
/*---------------------------------------------------------------------------*/
double BlockingRewardFunction::evaluate(std::shared_ptr<IState> AState) const {
	auto state = std::dynamic_pointer_cast<BlockingState>(AState);
	auto optimal = state->get_expected_optimal_score();
	auto memory = state->get_memory();
	auto reward = 0.0;

	reward = state->computeScore();

/*
	if (state->win()){
		return 1;
	}
	else if (state->lost())
		return -1;
		*/

	return reward;
/*
	if (state->win()){
		reward = state->get_expected_optimal_score();
	}
	else if (memory.size()==1){
		//root node
		reward = 0.0;
	}
	else{
		reward = *(memory.end()-1) - *(memory.end()-2);
	}
	return reward;*/
}
