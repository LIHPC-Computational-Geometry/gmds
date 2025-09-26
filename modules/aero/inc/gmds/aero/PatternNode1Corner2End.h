//
// Created by rochec on 28/04/23.
//

#ifndef GMDS_PATTERNNODE1CORNER2END_H
#define GMDS_PATTERNNODE1CORNER2END_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include <gmds/aero/AbstractPatternNode.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  PatternNode3Corner
 *  \brief
 */
class GMDSAero_API PatternNode1Corner2End: public AbstractPatternNode {
 public:
	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
    *  @param
	 */
	PatternNode1Corner2End(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
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

#endif     // GMDS_PATTERNNODE1CORNER2END_H
