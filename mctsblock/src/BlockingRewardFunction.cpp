/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/BlockingRewardFunction.h>
#include <gmds/mctsblock/BlockingState.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::mctsblock;
/*---------------------------------------------------------------------------*/
double BlockingRewardFunction::evaluate(std::shared_ptr<IState> AState) const {
	if(std::dynamic_pointer_cast<BlockingState>(AState)->win())
		return 1;
	else if(std::dynamic_pointer_cast<BlockingState>(AState)->lost())
		return -1;
	return 0;
}
