//
// Created by rochec on 28/10/2022.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_BEZIERSURFACE_H
#define GMDS_BEZIERSURFACE_H
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
     *  \brief Defines a Bezier surface in 3D.
 */
class GMDSMath_API BezierSurface {


 public:

	/*------------------------------------------------------------------------*/
	/** \brief Constructor of a bezier surface from an 2D array of ordered
	    * 	set of control points.
       *
       * \param APts the set of control points to define the curve
	 */
	BezierSurface(const Array2D<math::Point> & APts);

	/*------------------------------------------------------------------------*/
	/** \brief Returns the point located on the parametric surface with
       *
       * \param[in] u the first parameter in [0,1]
       * \param[in] v the second parameter in [0,1]
       *
       * \return The point located at (u,v) on the parametric surface
	 */
	math::Point operator()(const double& u, const double& v) const;

	/*------------------------------------------------------------------------*/

 private:
	/*------------------------------------------------------------------------*/
	/** \brief Returns the binomial coefficient C_k^n. This is the number of
	    *		part of k elements in a set of n elements
       *
       * \param[in] n
       * \param[in] k
       *
       * \return The binomial coefficient C_k^n
	 */
	double BinomialCoefficient(int n, int k) const;

	/*------------------------------------------------------------------------*/
	/** \brief Returns the value of the Bernstein's i-th polynomial of degree n
	    * 	for u in [0,1].
       *
       * \param[in] n degree of the Bernstein's polynomial
       * \param[in] i Bernstein's i-th polynomial of degree n
       * \param[in] u parameter in [0,1]
       *
       * \return Bernstein's i-th polynomial of degree n evaluated at u
	 */
	double BernsteinPolynomial(int n, int i, double u) const;

	/*------------------------------------------------------------------------*/
 private:
	/** ordered 2D array of control points of the curve */
	Array2D<math::Point> m_control_points;
	/** polynomial degrees of the surface */
	int m_degree_n;
	int m_degree_m;

};
/*----------------------------------------------------------------------------*/
}
  /*----------------------------------------------------------------------------*/
}

#endif     // GMDS_BEZIERSURFACE_H
