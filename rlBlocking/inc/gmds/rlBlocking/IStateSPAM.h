/*---------------------------------------------------------------------------*/
#ifndef MATCHING_ISTATE_H
#define MATCHING_ISTATE_H
/*---------------------------------------------------------------------------*/
#include <vector>
#include <memory>
#include <gmds/rlBlocking/IActionSPAM.h>
/*---------------------------------------------------------------------------*/
class IState {
public:
    virtual ~IState() = default;
    /**@brief Get the list of all the actions that can be applied on this state
     * @return a vector of actions to try
     */
    virtual std::vector<std::shared_ptr<IAction> > get_actions() const = 0;

    /**@brief Computes the state reach from the current one when applying @p AAction
     * @param AAction the action to apply
     * @return the state that is built from this one when applying @p AAction
     */
    virtual std::shared_ptr<IState> apply(std::shared_ptr<IAction> AAction) const = 0;

    virtual bool is_terminal() const = 0;
};


#endif //MATCHING_ISTATE_H
