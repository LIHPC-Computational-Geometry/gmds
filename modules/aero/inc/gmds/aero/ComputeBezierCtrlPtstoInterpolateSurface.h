//
// Created by rochec on 18/01/24.
//

#ifndef GMDS_COMPUTEBEZIERCTRLPTSTOINTERPOLATESURFACE_H
#define GMDS_COMPUTEBEZIERCTRLPTSTOINTERPOLATESURFACE_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/cadfac/FACManager.h>
#include <gmds/utils/Array.h>
#include <gmds/math/BezierSurface.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API ComputeBezierCtrlPtstoInterpolateSurface
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
         *  @param[in] AMeshHex the hex mesh we want to curve
         *
	 */
	ComputeBezierCtrlPtstoInterpolateSurface(cad::GeomSurface* ASurface,
	                                         Array2D<math::Point>* AInitCtrlPts);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief
	 *
	 * @return
	 */
	Array2D<math::Point> getCtrlPts();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Initialize the problem from the four nodes of the quad face.
	 * Will compute a set of points to interpolate on the surface.
	 *
	 * @return
	 */
	void InitInterpolatePoints();
	/*-------------------------------------------------------------------*/
	/** @brief Compute control points in order to interpolate some points
	 * on the geometry/surface.
	 *
	 * @return
	 */
	void ComputeControlPoints();
	/*-------------------------------------------------------------------*/

 private:
	/** Geometric surface to interpolate */
	cad::GeomSurface* m_Surface;
	/** Degree m */
	int m_degree_m;
	/** Degree n */
	int m_degree_n;
	/** Init Ctrl Pts */
	Array2D<math::Point>* m_InitCtrlPts;
	/** Ctrl Pts */
	Array2D<math::Point>* m_CtrlPts;
	/** Points to interpolate */
	Array2D<math::Point>* m_InterpolatePts;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_COMPUTEBEZIERCTRLPTSTOINTERPOLATESURFACE_H
