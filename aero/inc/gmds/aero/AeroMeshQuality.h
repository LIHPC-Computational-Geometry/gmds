//
// Created by rochec on 14/04/2022.
//

#ifndef GMDS_AEROMESHQUALITY_H
#define GMDS_AEROMESHQUALITY_H

/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {
/*----------------------------------------------------------------------------*/
/** \class Utils
 *  \brief
 **/
class LIB_GMDS_AERO_API AeroMeshQuality {
 public:
	/*------------------------------------------------------------------------*/
	/** \brief  The angles side by side the edge (n0_id, n1_id)
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double minlenghtedge(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  The sum of the angle between the vectors p0p1 and p0P2, and the
	 * vectors p0p1 and p0p3
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the angle
	 */
	static double AngleOuverture(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Aspect Ratio
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double AspectRatioQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Internal Angle Deviation
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double InternalAngleDeviationQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Equi Angle Skewness
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the min ratio between two opposite edges (between 0 and 1)
	 */
	static double EquiAngleSkewnessQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Condition for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* Définie entre 1 et double_max, une valeur entre 1 et 4 est acceptable.
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the condition of the quad
	 */
	static double ConditionQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Edge ratio for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* Définie entre 1 et double_max, une valeur entre 1 et 1.3 est acceptable.
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the edge ratio of the quad
	 */
	static double EdgeRatioQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Jacobian for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* Définie entre 0 et double_max, une valeur entre 0 et double_max est acceptable.
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the jacobian of the quad
	 */
	static double JacobianQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Scaled Jacobian for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* Définie entre -1 et 1, une valeur entre 0.3 et 1 est acceptable.
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the scaled jacobian of the quad
	 */
	static double ScaledJacobianQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Shape for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* Définie entre 0 et 1, une valeur entre 0.3 et 1 est acceptable.
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the shape of the quad
	 */
	static double ShapeQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Skew for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* Définie entre 0 et 1, une valeur entre 0 et 0.5 est acceptable.
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
         *
         * \return  the skew of the quad
	 */
	static double SkewQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Stretch for QUAD (cf The Verdict Geometric Quality Library from SANDIA)
	 		* Définie entre 0 et 1, une valeur entre 0.25 et 1 est acceptable.
         *
         * \param[in] p0 first point
         * \param[in] p1 second point
         * \param[in] p2 third point
         * \param[in] p3 fourth point
         *
         * \return  the stretch of the quad
	 */
	static double StretchQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3);
	/*------------------------------------------------------------------------*/
	/** \brief  Scaled Jacobian for HEX (cf The Verdict Geometric Quality Library from SANDIA)
	 		* This metric is the minimum determinant of the Jacobian matrix evaluated at each
	 		* corner and the center of the element, divided by the corresponding edge lengths.
	 		* Defined between -1 et 1, a value between 0.5 and 1 is acceptable.
         *
         * \param[in] p0 first point
			* \param[in] p1 second point
			* \param[in] p2 third point
			* \param[in] p3 fourth point
			* \param[in] p4 first point
			* \param[in] p5 second point
			* \param[in] p6 third point
			* \param[in] p7 fourth point
         *
         * \return the scaled jacobian of the hex
	 */
	static double ScaledJacobianHEX(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3,
	                                const math::Point& p4, const math::Point& p5, const math::Point& p6, const math::Point& p7);
	/*------------------------------------------------------------------------*/


};
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROMESHQUALITY_H
