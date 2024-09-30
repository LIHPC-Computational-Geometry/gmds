//
// Created by rochec on 13/02/24.
//

#ifndef GMDS_COMPUTEBEZIERCURVECTRLPTSTOINTERPOLATECURVE_H
#define GMDS_COMPUTEBEZIERCURVECTRLPTSTOINTERPOLATECURVE_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/cadfac/FACManager.h>
#include <gmds/utils/Array.h>
#include <gmds/math/BezierSurface.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API ComputeBezierCurveCtrlPtstoInterpolateCurve
{

 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param[in] ACurve the geometric curve we want to fit
         *  @param[in] AInitCtrlPts the init control points of the Bézier Curve
         *
	 */
	ComputeBezierCurveCtrlPtstoInterpolateCurve(cad::GeomCurve* ACurve,
	                                            std::vector<math::Point>* AInitCtrlPts);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief Return the vector of the control points of the Bézier curve.
	 *
	 * @return
	 */
	std::vector<math::Point> getCtrlPts();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Initialize the problem from the two bounding points of the
	 * curve.
	 * Will compute a set of points to interpolate the curve.
	 *
	 * @return
	 */
	void InitInterpolatePoints();
	/*-------------------------------------------------------------------*/
	/** @brief Compute control points in order to interpolate some points
	 * on the geometric curve.
	 *
	 * @return
	 */
	void ComputeControlPoints();
	/*-------------------------------------------------------------------*/

 private:
	/** Geometric surface to interpolate */
	cad::GeomCurve* m_Curve;
	/** Degree n */
	int m_degree_n;
	/** Init Ctrl Pts */
	std::vector<math::Point>* m_InitCtrlPts;
	/** Ctrl Pts */
	std::vector<math::Point>* m_CtrlPts;
	/** Points to interpolate */
	std::vector<math::Point>* m_InterpolatePts;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_COMPUTEBEZIERCURVECTRLPTSTOINTERPOLATECURVE_H
