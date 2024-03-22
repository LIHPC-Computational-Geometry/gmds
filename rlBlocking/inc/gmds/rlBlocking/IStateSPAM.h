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

	 /**@brief Indicates if the concrete state is terminal (win or lost in classical games)
     * @return true if it is terminal, false otherwise
	  */
    virtual bool is_terminal() const = 0;

	 /**@provide a function to write a file that stores the state knowing a file name prefix @p AFileName,
     * the stage index @p AStageIndex, the id, and the depth of the node that knows this state in the
     * MCTS tree. By default, you can provide an empty implementation that does nothing.
     *
     * @param[in] AFileName  A file name
     * @param[in] AStageIndex Stage index
     * @param[in] ANodeId    Id of the node that contains this state
     * @param[in] ADepth     Depth of the node that contains this state
	  */
	 virtual void write(const std::string& AFileName,
	                    const int AStageIndex,
	                    const int ANodeId,
	                    const int ADepth) const = 0;
};


#endif //MATCHING_ISTATE_H
