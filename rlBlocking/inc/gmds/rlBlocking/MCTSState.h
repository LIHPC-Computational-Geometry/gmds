//
// Created by bourmaudp on 02/12/22.
//
/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSSTATE_H
#define GMDS_MCTSSTATE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSMove.h>
#include <iostream>
/*----------------------------------------------------------------------------------------*/
#include <queue>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSState
 *  @brief  Class that provides the interface to be implemented for performing the
 *  MCST algorithm
 */
class LIB_GMDS_RLBLOCKING_API MCTSState {
 public:
		/*--------------------------------------------------------------------*/
		/** @enum  Status code for rollout execution
		 */
		typedef enum {
			WIN,
			LOSE,
			DRAW
		} ROLLOUT_STATUS;
	   /*------------------------------------------------------------------------*/
	/** @brief  Destructor
	 */
	virtual ~MCTSState() = default;
	/*------------------------------------------------------------------------*/
	/** @brief  Gives the set of actions that can be tried from the current state
	 */
	virtual std::deque<MCTSMove *> *actions_to_try() const = 0;
	/*------------------------------------------------------------------------*/
	/** @brief  Performs the @p AMove to change of states
	 * @param[in] AMove the movement to apply to get to a new state
	 */
	virtual MCTSState *next_state(const MCTSMove *AMove) const = 0;
	/*------------------------------------------------------------------------*/
	/** @brief Rollout from this state (random simulation)
	 *  @return the rollout status
	 */
	virtual double state_rollout() const = 0;
	/*------------------------------------------------------------------------*/
	/** @brief check the result of a terminal state
	 *  @return the value of the result: Win, Lose, Draw
	 */
	virtual ROLLOUT_STATUS result_terminal() const = 0;
	/*------------------------------------------------------------------------*/
	/** @brief  Indicate if we have a terminal state (win=true, fail=false)
	 * @return true if we have a leaf (in the sense of a traditional tree)
	 */
	virtual bool is_terminal() const = 0;
	/*------------------------------------------------------------------------*/
	/** @brief  Indicate if we have a terminal state (win=true, fail=false)
	 * @return true if we have a leaf (in the sense of a traditional tree)
	 */
	virtual double get_quality() const = 0;

	virtual void print() const {
		std::cout << "Printing not implemented" << std::endl;
	}
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSSTATE_H
/*----------------------------------------------------------------------------------------*/
