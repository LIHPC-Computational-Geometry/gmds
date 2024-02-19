/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSMOVE_POLYCUBE_H
#define GMDS_MCTSMOVE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSMove.h>
#include <gmds/utils/CommonTypes.h>
#include <deque>
#include <random>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSMove
 *  @brief  Structure that provides ....
 */
struct LIB_GMDS_RLBLOCKING_API MCTSMovePolycube: public MCTSMove {
	/*------------------------------------------------------------------------*/
	/** @brief  Destructor
	 */
	~MCTSMovePolycube();
	/*------------------------------------------------------------------------*/
	TCellID m_AIdEdge;
	TCellID m_AIdBlock;
	double m_AParamCut;
	/** @brief if typeMove=2: delete block, typeMove=1 cut block, typeMove=3 updateClassification,
	 */
	unsigned int m_typeMove;

	/** @brief  Overloaded ==
	 */
	MCTSMovePolycube(TCellID AIdEdge = -1,TCellID AIdBlock = -1 , double AParamCut = 0,unsigned int ATypeMove = -1);
	bool operator==(const MCTSMove& AOther) const;
	void print() const;



};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSMOVE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
