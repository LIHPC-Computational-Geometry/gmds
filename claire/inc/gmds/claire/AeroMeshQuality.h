//
// Created by rochec on 14/04/2022.
//

#ifndef GMDS_AEROMESHQUALITY_H
#define GMDS_AEROMESHQUALITY_H

/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {
/*----------------------------------------------------------------------------*/
/** \class Utils
 *  \brief
 **/
class LIB_GMDS_CLAIRE_API AeroMeshQuality {
 public:
	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance between two nodes given the ids
         *
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double oppositeedgeslenghtratio(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);

	/*------------------------------------------------------------------------*/
	/** \brief  The angles side by side the edge (n0_id, n1_id)
         *
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double angleouverture(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  The angles side by side the edge (n0_id, n1_id)
         *
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double minlenghtedge(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Aspect Ratio
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double AspectRatioQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Internal Angle Deviation
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double InternalAngleDeviationQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Equi Angle Skewness
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double EquiAngleSkewnessQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/



};
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROMESHQUALITY_H
