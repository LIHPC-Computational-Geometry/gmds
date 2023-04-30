//
// Created by rochec on 28/04/23.
//

#ifndef GMDS_PATTERNNODE3END_H
#define GMDS_PATTERNNODE3END_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractPatternNode.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  PatternNode3Corner
 *  \brief
 */
class LIB_GMDS_CLAIRE_API PatternNode3End: public AbstractPatternNode {
 public:
	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
    *  @param
	 */
	PatternNode3End(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
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

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_PATTERNNODE3END_H
