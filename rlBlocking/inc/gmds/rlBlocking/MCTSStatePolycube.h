//
// Created by bourmaudp on 02/12/22.
//
/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSSTATE_POLYCUBE_H
#define GMDS_MCTSSTATE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSState.h>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSState
 *  @brief  Class that provides the interface to be implemented for performing the
 *  MCST algorithm
 */
class LIB_GMDS_RLBLOCKING_API MCTSStatePolycube: public MCTSState{
 public:
	   /*------------------------------------------------------------------------*/
	/** @brief  Destructor
	 */
	virtual ~MCTSStatePolycube();
	/*------------------------------------------------------------------------*/
	/** @brief  Gives the set of actions that can be tried from the current state
	 */
	virtual std::queue<MCTSMove *> *actions_to_try() const ;
	/*------------------------------------------------------------------------*/
	/** @brief  Performs the @p AMove to change of states
	 * @param[in] AMove the movement to apply to get to a new state
	 */
	virtual MCTSState *next_state(const MCTSMove *AMove) const;
	/*------------------------------------------------------------------------*/
	/** @brief Rollout from this state (random simulation)
	 *  @return the rollout status
	 */
	virtual ROLLOUT_STATUS rollout() const;
	/*------------------------------------------------------------------------*/
	/** @brief  Indicate if we have a terminal state (win=true, fail=false)
	 * @return true if we have a leaf (in the sense of a traditional tree)
	 */
	virtual bool is_terminal() const;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSSTATE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
