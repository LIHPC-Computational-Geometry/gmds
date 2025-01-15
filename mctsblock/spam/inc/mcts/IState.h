/*---------------------------------------------------------------------------*/
#ifndef MATCHING_ISTATE_H
#define MATCHING_ISTATE_H
/*---------------------------------------------------------------------------*/
#include <vector>
#include <memory>
#include <mcts/IAction.h>
/*---------------------------------------------------------------------------*/
class IState {
public:
    virtual ~IState() = default;
    /**@brief Get the list of all the actions that can be applied on this state
     * @return a vector of actions to try
     */
    virtual std::vector<std::shared_ptr<IAction> > get_actions_selection()  = 0;
	/**@brief Get the limited list of  the actions that can be applied on this state.
	 *Impossible to remove blocks in the geometric volume. Limit the possible cuts.
	 * @return a vector of actions to try
	 */
    virtual std::vector<std::shared_ptr<IAction> > get_actions_simulation()  = 0;

    /**@brief Indicates if the concrete state is terminal (win or lost in classical games)
     * @return true if it is terminal, false otherwise
     */
    virtual bool is_terminal()  = 0;

		  /**@brief Indicates if the concrete state is a winning state
			 * @return true if it is, false otherwise
			*/
		  virtual bool win()  = 0;
		  /**@brief Indicates if the concrete state is a losing state
			 * @return true if it is, false otherwise
			*/
		  virtual bool lost()  = 0;
	     /**@brief Indicates if the concrete state is a draw state
      * @return true if it is, false otherwise
	      */
	     virtual bool draw()  = 0;

    /**@provide a function to write a file that stores the state knowing a file name prefix @p AFileName,
     * the stage index @p AStageIndex, the id, and the depth of the node that knows this state in the
     * MCTS tree. By default, you can provide an empty implementation that does nothing.
     * The returned value is a string that will be stored in a "data" field of the json tree output. By default
     * return an empty string if you don't want to store something. Oterwise, please respect json syntax.
     *
     * @param[in] AFileName  A file name
     * @param[in] AStageIndex Stage index
     * @param[in] ANodeId    Id of the node that contains this state
     * @param[in] ADepth     Depth of the node that contains this state
     *
     * @return A string to store in the data field of each node
     */
    virtual std::string write(const std::string& AFileName,
                       const int AStageIndex,
                       const int ANodeId,
                       const int ADepth) const = 0;
};


#endif //MATCHING_ISTATE_H
