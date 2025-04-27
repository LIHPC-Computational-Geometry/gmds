/*---------------------------------------------------------------------------*/
#ifndef MATCHING_IREWARDFUNCTION_H
#define MATCHING_IREWARDFUNCTION_H
/*---------------------------------------------------------------------------*/
#include <memory>
#include <mcts/IAction.h>
#include <mcts/IState.h>
/*---------------------------------------------------------------------------*/
struct IRewardFunction {
    virtual ~IRewardFunction() = default;

    /**@brief this function must give the expected reward for
     * state @p AState
     * @param AState
     * @return a double value
     */
    virtual double evaluate(std::shared_ptr<IState> AState) const = 0;
};
/*---------------------------------------------------------------------------*/
#endif //MATCHING_IREWARDFUNCTION_H
/*---------------------------------------------------------------------------*/
