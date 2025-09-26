/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_BLOCKING_REWARD_FUNCTION_H
#define GMDS_MCTS_BLOCKING_REWARD_FUNCTION_H
/*----------------------------------------------------------------------------*/
#include <GMDSMctsBlock_export.h>
#include "mcts/IRewardFunction.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace mctsblock {
/*----------------------------------------------------------------------------*/
class GMDSMctsBlock_API BlockingRewardFunction: public IRewardFunction{

	double evaluate(std::shared_ptr<IState> AState) const override;

};
/*----------------------------------------------------------------------------*/
}     // namespace blocking
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MCTS_BLOCKING_REWARD_FUNCTION_H
