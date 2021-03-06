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
	/** \brief  Condition for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* D??finie entre 1 et double_max, une valeur entre 1 et 4 est acceptable.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the condition of the quad
	 */
	static double ConditionQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Edge ratio for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* D??finie entre 1 et double_max, une valeur entre 1 et 1.3 est acceptable.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the edge ratio of the quad
	 */
	static double EdgeRatioQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Jacobian for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* D??finie entre 0 et double_max, une valeur entre 0 et double_max est acceptable.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the jacobian of the quad
	 */
	static double JacobianQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Scaled Jacobian for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* D??finie entre -1 et 1, une valeur entre 0.3 et 1 est acceptable.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the scaled jacobian of the quad
	 */
	static double ScaledJacobianQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Shape for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* D??finie entre 0 et 1, une valeur entre 0.3 et 1 est acceptable.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the shape of the quad
	 */
	static double ShapeQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Skew for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* D??finie entre 0 et 1, une valeur entre 0 et 0.5 est acceptable.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the skew of the quad
	 */
	static double SkewQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/
	/** \brief  Stretch for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* D??finie entre 0 et 1, une valeur entre 0.25 et 1 est acceptable.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id third node id
         * \param[in] n3_id fourth node id
         *
         * \return  the stretch of the quad
	 */
	static double StretchQUAD(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*------------------------------------------------------------------------*/


};
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROMESHQUALITY_H
