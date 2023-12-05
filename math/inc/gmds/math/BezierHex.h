//
// Created by rochec on 05/12/23.
//

#ifndef GMDS_BEZIERHEX_H
#define GMDS_BEZIERHEX_H
/*----------------------------------------------------------------------------*/
#include <cmath>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Array.h>
#include <gmds/math/Point.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*--------------------------------------------------------------------------*/
namespace math {
/*------------------------------------------------------------------------*/
/** \class BezierSurface
     *  \brief Defines a Bezier hex in 3D.
 */
class GMDSMath_API BezierHex {


 public:

	/*------------------------------------------------------------------------*/
	/** \brief Constructor of a bezier hex from an 3D array of ordered
	    * 	set of control points.
       *
       * \param[in] APts the set of control points to define the curve
	 */
	BezierHex(const Array3D<math::Point> & APts);

	/*------------------------------------------------------------------------*/
	/** \brief Returns the point located on the parametric surface with:
       *
       * \param[in] Au the first parameter in [0,1]
       * \param[in] Av the second parameter in [0,1]
       * \param[in] Aw the third parameter in [0,1]
       *
       * \return The point located at (Au,Av,Aw) on the parametric surface.
	 */
	math::Point operator()(const double& Au, const double& Av, const double& Aw) const;

	/*------------------------------------------------------------------------*/
	/** \brief Returns a set of point that discretize the volume in
	    * 	ANb_l x ANb_m x ANb_n segments in the parametric space.
       *
       * \param[in] ANb_l the number of segments in the first direction
       * \param[in] ANb_m the number of segments in the second direction
       * \param[in] ANb_n the number of segments in the third direction
       *
       * \return The set of point discretizing the volume (*this).
	 */
	Array3D<math::Point> getDiscretization(const int ANb_l, const int ANb_m, const int ANb_n) const;

 private:
	/*------------------------------------------------------------------------*/
	/** \brief Returns the binomial coefficient C_Ak^An. This is the number of
	    *		part of Ak elements in a set of An elements.
       *
       * \param[in] An the number of elements in the set
       * \param[in] Ak
       *
       * \return The binomial coefficient C_Ak^An.
	 */
	double BinomialCoefficient(const int An, const int Ak) const;

	/*------------------------------------------------------------------------*/
	/** \brief Returns the value of the Bernstein's Ai-th polynomial of degree An
	    * 	for Au in [0,1].
       *
       * \param[in] An degree of the Bernstein's polynomial
       * \param[in] Ai Bernstein's i-th polynomial of degree n
       * \param[in] Au parameter in [0,1]
       *
       * \return Bernstein's Ai-th polynomial of degree An evaluated at Au
	 */
	double BernsteinPolynomial(const int An, const int Ai, const double Au) const;

 private:
	/** ordered 3D array of control points of the volume */
	Array3D<math::Point> m_control_points;
	/** polynomial degree of the volume in the first direction */
	int m_degree_l;
	/** polynomial degree of the volume in the second direction */
	int m_degree_m;
	/** polynomial degree of the surface in the third direction */
	int m_degree_n;

};
/*----------------------------------------------------------------------------*/
}
  /*----------------------------------------------------------------------------*/
}
#endif     // GMDS_BEZIERHEX_H
