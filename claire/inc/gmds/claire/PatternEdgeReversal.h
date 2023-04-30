//
// Created by rochec on 28/04/23.
//

#ifndef GMDS_PATTERNEDGEREVERSAL_H
#define GMDS_PATTERNEDGEREVERSAL_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractPatternEdge.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  PatternEdgeReversal
 *  \brief
 */
class LIB_GMDS_CLAIRE_API PatternEdgeReversal: public AbstractPatternEdge {
 public:
	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
    *  @param
	 */
	PatternEdgeReversal(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
	               Mesh *AMeshT, FastLocalize *Afl,
	               double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField);
	/*-------------------------------------------------------------------*/
 protected:
	/*-------------------------------------------------------------------*/
	/** \brief
	 * @param
	 * @return
	 */
	void computeNewHex() override;
	/*-------------------------------------------------------------------*/
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_PATTERNEDGEREVERSAL_H
