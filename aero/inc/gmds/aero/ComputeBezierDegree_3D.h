//
// Created by rochec on 17/01/24.
//

#ifndef GMDS_COMPUTEBEZIERDEGREE_3D_H
#define GMDS_COMPUTEBEZIERDEGREE_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/cadfac/FACManager.h>
#include <gmds/utils/Array.h>
#include <gmds/math/BezierSurface.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API ComputeBezierDegree_3D
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
         *  @param[in] AManager geometric manager
         *  @param[in] ALinker_HG linker of the hex mesh to the geom
         *  @param[in] Afaces_to_fit list of the faces we want to fit to the geom
         *  @param[in] Amax_error error tolerance
         *  @param[in] Amax_degree max degree allowed
         *
	 */
	ComputeBezierDegree_3D(Mesh *AMeshHex, cad::FACManager* AManager, cad::GeomMeshLinker* ALinker_HG, std::vector<TCellID>* Afaces_to_fit, double Amax_error, int Amax_degree);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief
	 *
	 * @return the final degree computed to respect the max error.
	 */
	int getDegree();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the control points of the Bezier surface of degree
	 * degree x degree in order to interpolate the geometric surface surf.
	 *
	 * @param[in] f Face
	 * @param[in] degree Degree of the Bezier surface
	 * @param[in] surf Geometric surface to fit
	 *
	 * @return the length of the bezier curve
	 */
	Array2D<math::Point> computeBezierCtrlPtstoInterpolateSurface(Face f, int degree, cad::GeomSurface* surf);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the error between the Bezier surface Abs and the
	 * geometric surface surf. The error computed is the maximal norm between
	 * a point on the Bezier surface and its projection onto the geometric
	 * surface surf. The error is computed on a sample of 20x20 points taken
	 * uniformly in the parametric space of the Bezier Surface.
	 *
	 * @param[in] Abs Bezier Surface
	 *
	 * @return the length of the bezier curve
	 */
	double computeErrorBtwBezierSurfaceandGeomSurface(math::BezierSurface Abs, cad::GeomSurface* surf);
	/*-------------------------------------------------------------------*/

 private:
	/** hex mesh to curve */
	Mesh *m_meshHex;
	/** Manager */
	cad::FACManager* m_manager;
	/** Linker Hex mesh to geometry */
	cad::GeomMeshLinker* m_linker_HG;
	/** List of the faces to fit to the geom */
	std::vector<TCellID>* m_faces_to_fit;
	/** Error tolerance */
	double m_max_error;
	/** Max degree allowed */
	int m_max_degree;
	/** final degree computed */
	int m_final_degree;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_COMPUTEBEZIERDEGREE_3D_H
