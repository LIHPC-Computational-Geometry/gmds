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
	/** \brief Compute the normal vector to the faces around the node @p n_id
	 * of the front, but only on the same side than the face @p f_id. We
	 * consider here that around each node of the edge m_e_id, there are two
	 * sides, separeted by the corner edges path.
	 *
	 * @param[in] n_id id of one node of the corner edge @p m_e_id
	 * @param[in] f_id id of one face adjacent to the corner edge @p m_e_id to
	 * 			define the side
	 *
	 * @return A 3D vector corresponding to the mean of the normals of the faces
	 * 		around @p n_id of the side @p f_id
	 */
	math::Vector3d computeNormaltoFacesAroundNodeSideFace(TCellID n_id, TCellID f_id);
	/*-------------------------------------------------------------------*/
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_PATTERNEDGECORNER_H
