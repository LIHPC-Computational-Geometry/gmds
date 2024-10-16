//
// Created by rochec on 27/04/23.
//

#ifndef GMDS_PATTERNNODE3CORNER_H
#define GMDS_PATTERNNODE3CORNER_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/aero/AbstractPatternNode.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  PatternNode3Corner
 *  \brief
 */
class LIB_GMDS_AERO_API PatternNode3Corner: public AbstractPatternNode {
 public:
	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
    *  @param
	 */
	PatternNode3Corner(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
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

#endif     // GMDS_PATTERNNODE3CORNER_H
