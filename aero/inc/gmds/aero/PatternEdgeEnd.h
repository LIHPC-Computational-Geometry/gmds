//
// Created by rochec on 28/04/23.
//

#ifndef GMDS_PATTERNEDGEEND_H
#define GMDS_PATTERNEDGEEND_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/aero/AbstractPatternEdge.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  PatternEdgeEnd
 *  \brief
 */
class LIB_GMDS_AERO_API PatternEdgeEnd: public AbstractPatternEdge {
 public:
	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
    *  @param
	 */
	PatternEdgeEnd(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
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

#endif     // GMDS_PATTERNEDGEEND_H
