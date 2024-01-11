/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSMOVE_POLYCUBE_H
#define GMDS_MCTSMOVE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSMove.h>
#include <gmds/utils/CommonTypes.h>
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
	/** @brief if typeMove=0: delete block, typeMove=1 cut block
	 */
	bool m_typeMove;

	/** @brief  Overloaded ==
	 */
	MCTSMovePolycube(TCellID AIdEdge,TCellID AIdBlock, double AParamCut,bool ATypeMove);
	bool operator==(const MCTSMove& AOther) const;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSMOVE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
