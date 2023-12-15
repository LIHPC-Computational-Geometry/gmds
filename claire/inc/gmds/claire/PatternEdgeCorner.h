//
// Created by rochec on 28/04/23.
//

#ifndef GMDS_PATTERNEDGECORNER_H
#define GMDS_PATTERNEDGECORNER_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractPatternEdge.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  PatternEdgeCorner
 *  \brief
 */
class LIB_GMDS_CLAIRE_API PatternEdgeCorner: public AbstractPatternEdge {
 public:
	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
    *  @param
	 */
	PatternEdgeCorner(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
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
	/** \brief
	 * @param
	 * @return
	 */
	math::Vector3d computeNormaltoFacesAroundNodeSideFace(TCellID n_id, TCellID f_id);
	/*-------------------------------------------------------------------*/
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_PATTERNEDGECORNER_H
