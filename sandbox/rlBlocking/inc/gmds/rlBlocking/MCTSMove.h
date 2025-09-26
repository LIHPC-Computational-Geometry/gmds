//
// Created by bourmaudp on 02/12/22.
//
/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSMOVE_H
#define GMDS_MCTSMOVE_H
/*----------------------------------------------------------------------------------------*/
#include "GMDSRlBlocking_export.h"
#include <string>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSMove
 *  @brief  Structure that provides ....
 */
struct GMDSRlBlocking_API MCTSMove {
	/*------------------------------------------------------------------------*/
	/** @brief  Destructor
	 */
	virtual ~MCTSMove() = default;
	/*------------------------------------------------------------------------*/
	/** @brief  Overloaded ==
	 */
	virtual bool operator==(const MCTSMove& AOther) const = 0;
	virtual std::string sprint() const { return "Not implemented"; }
	virtual void print() const =0;   // and optionally this
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSMOVE_H
/*----------------------------------------------------------------------------------------*/
